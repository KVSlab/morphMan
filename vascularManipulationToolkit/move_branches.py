from common import *
from argparse import ArgumentParser
from os import path, listdir
from subprocess import STDOUT, check_output
from time import time
from copy import deepcopy

from IPython import embed
import sys
import math

# Local import
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
import ToolRepairSTL

def read_command_line():
    """
    Read arguments from commandline
    """
    # TODO: Change to --feature --no-feature for bool agruments
    # Ex. 
    # feature_parser = parser.add_mutually_exclusive_group(required=False)
    # feature_parser.add_argument('--feature', dest='feature', action='store_true')
    # feature_parser.add_argument('--no-feature', dest='feature', action='store_false')
    # parser.set_defaults(feature=True)

    # Also add choise when that is approporiate

    # Could also make bool be int and add choise 0, 1 to be true false...
    parser = ArgumentParser()

    parser.add_argument('-d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('-c', '--case', type=str, default=None, help="Choose case")
    parser.add_argument('--s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram (surface) should be" + \
                        "smoothed before it is manipulated", metavar="smooth")
    parser.add_argument('--a', '--angle', type=float, default=10,
                        help="Each daughter branch is rotated an angle a in the" + \
                        " bifurcation plane. a should be expressed in radians as" + \
                        " any math expression from the" + \
                        " math module in python", metavar="rotation_angle")
    parser.add_argument('--smooth_factor', type=float, default=0.25,
                         help="If smooth option is true then each voronoi point" + \
                         " that has a radius less then MISR*(1-smooth_factor) at" + \
                         " the closest centerline point is removes",
                         metavar="smoothening_factor")
    parser.add_argument("--leave1", type=bool, default=False,
                        help="Leave one branch untuched")
    parser.add_argument("--leave2", type=bool, default=False,
                        help="Leave one branch untuched")
    parser.add_argument("--bif", type=bool, default=False,
                        help="interpolate bif as well")
    parser.add_argument("--lower", type=bool, default=False,
                        help="Make a fourth line to interpolate along that" + \
                             " is lower than the other bif line.")
    parser.add_argument("--cylinder_factor", type=float, default=7.0,
                        help="Factor for choosing the smaller cylinder")
    parser.add_argument("--version", type=bool, default=True, help="Type of" + \
                        "interpolation")
    parser.add_argument("--aneurysm", type=bool, default=False,
                        help="Determines if there is an aneurysm or not")
    parser.add_argument("--anu_num", type=int, default=0,
                        help="If multiple aneurysms, choise one")

    parser.add_argument("--step", type=float, default=0.1,
                        help="Resampling step used to resample centerlines")

    args = parser.parse_args()
    ang_ = 0 if args.a == 0 else math.pi / args.a

    return args.s, ang_, args.smooth_factor, args.leave1, args.leave2, \
           args.bif, args.dir_path, args.case, args.lower, \
           args.cylinder_factor, args.version, args.aneurysm, args.anu_num, args.step



def get_points(data, key, R, m, rotated=True, bif=False):
    """
    Finds spesific points around the bifurcation, based on the
    key argument. Points can before or after rotation. 

    Args:
        data (dict): Contains information about points and IDs of branches and bifurcation.
        key (str): Type of points to extract.
        R (ndarray): Matrix containing unit vectors in the rotated coordinate system.
        m (dict): Cointains rotation matrices for each daughter branch.
        rotated (bool): Gets rotated points if True.
        bif (true): Gets only bifurcation points if True.
    
    Returns:
        points (vtkPoints): Points as VTK objects.
    Returns:
        div_points_bif (ndarray): Points as numpy objects. 
    """
    div_points = np.asarray([data["bif"][key], data[0][key], data[1][key]])

    # Origo of the bifurcation
    O_key = "div_point"
    O = np.asarray([data["bif"][O_key], data[0][O_key], data[1][O_key]])
    O = np.sum(np.asarray(O),axis=0)/3.

    if rotated:
        R_inv = np.linalg.inv(R)
        for i in range(len(div_points)):
            m_ = m[i] if i > 0 else np.eye(3)
            div_points[i] = np.dot(np.dot(np.dot(div_points[i] - O, R), m_), R_inv) + O

    # Insert landmarking points into VTK objects
    points = vtk.vtkPoints()
    div_points_bif = div_points[bif:]
    for point in div_points_bif:
        points.InsertNextPoint(point)
    
    return points, div_points_bif


def rotate_voronoi(clipped_voronoi, patch_cl, div_points, m, R):
    """
    Perform rotation of the voronoi diagram representing the 
    daughter branches. Rotate along the bifurcation plane
    spanned by two vectors, preserving the angle with
    the rest of the vasculature. Rotation is performed
    using a standard rotational matrix m. 

    Args:
        clipped_voronoi (vtkPolyData): Clipped voronoi diagram.
        patch_cl (vtkPolyData): Clipped centerline.
        div_points (ndarray): Contains bifurcation landmarking points. 
        R (ndarray): Matrix containing unit vectors in the rotated coordinate system.
        m (dict): Cointains rotation matrices for each daughter branch.
    Returns:
        maskedVoronoi (vtkPolyData): Rotated voronoi diagram. 
    """
    numberOfPoints = clipped_voronoi.GetNumberOfPoints()
    distance = vtk.vtkMath.Distance2BetweenPoints
    I = np.eye(3)
    R_inv = np.linalg.inv(R)

    locator = []
    cellLine = []
    not_rotate = [0]
    for i in range(patch_cl.GetNumberOfCells()):
        cellLine.append(extract_single_line(patch_cl, i))
        tmp_locator = get_locator(cellLine[-1])
        locator.append(tmp_locator)

    for i in range(1, patch_cl.GetNumberOfCells()):
        pnt = cellLine[i].GetPoints().GetPoint(0)
        new = cellLine[0].GetPoints().GetPoint(locator[0].FindClosestPoint(pnt))
        dist = math.sqrt(distance(pnt, new)) < divergingRatioToSpacingTolerance
        if dist:
            not_rotate.append(i)

    def check_rotate(point):
        dist = []
        for i in range(len(locator)):
            tmp = locator[i].FindClosestPoint(point)
            tmp = cellLine[i].GetPoints().GetPoint(tmp)
            dist.append(math.sqrt(distance(tmp, point)))

        if dist.index(min(dist)) not in not_rotate:
            pnt = cellLine[dist.index(min(dist))].GetPoints().GetPoint(0)
            if math.sqrt(distance(pnt, div_points[1])) >  \
                    math.sqrt(distance(pnt, div_points[2])):
                m_ = m[2]
                div = div_points[2]
            else:
                m_ = m[1]
                div = div_points[1]
            return m_, div
        else:
            return I, np.array([0, 0, 0])
    
    maskedVoronoi = vtk.vtkPolyData()
    maskedPoints = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPoints) 

    # Iterate through voronoi diagram
    for i in range(numberOfPoints):
        point = [0.0, 0.0, 0.0]
        clipped_voronoi.GetPoint(i, point)

        pointRadius = clipped_voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        M, O = check_rotate(point)
        tmp = np.dot(np.dot(np.dot(np.asarray(point) - O, R), M), R_inv) + O
        maskedPoints.InsertNextPoint(tmp)
        radiusArray.SetTuple1(i, pointRadius)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i)

    maskedVoronoi.SetPoints(maskedPoints)
    maskedVoronoi.SetVerts(cellArray)
    maskedVoronoi.GetPointData().AddArray(radiusArray)

    return maskedVoronoi


def rotate_cl(patch_cl, div_points, rotation_matrix, R):
    """
    Perform rotation of the centerline representing the 
    daughter branches. Rotate along the bifurcation plane
    spanned by two vectors, preserving the angle with
    the rest of the vasculature. Rotation is performed
    using a standard rotational matrix. 

    Args:
        patch_cl (vtkPolyData): Clipped centerline representing two daughter branches.
        div_points (ndarray): Contains bifurcation landmarking points. 
        rotation_matrix (dict): Cointains rotation matrices for each daughter branch.
        R (ndarray): Matrix containing unit vectors in the rotated coordinate system.
    Returns:
        centerlin (vtkPolyData): Rotated centerline. 
    """
    distance = vtk.vtkMath.Distance2BetweenPoints
    I = np.eye(3)
    R_inv = np.linalg.inv(R)
    
    numberOfPoints = patch_cl.GetNumberOfPoints()

    centerline = vtk.vtkPolyData()
    centerlinePoints = vtk.vtkPoints()
    centerlineCellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPoints)

    line0 = extract_single_line(patch_cl, 0)
    locator0 = get_locator(line0)

    # Iterate through points along the centerline 
    count = 0
    for i in range(patch_cl.GetNumberOfCells()):
        cell = extract_single_line(patch_cl, i)
        centerlineCellArray.InsertNextCell(cell.GetNumberOfPoints())

        start = cell.GetPoint(0)
        dist = line0.GetPoint(locator0.FindClosestPoint(start))
        test = math.sqrt(distance(start, dist)) > divergingRatioToSpacingTolerance

        if test or len(div_points) == 2:
            locator = get_locator(cell)

            pnt1 = cell.GetPoint(locator.FindClosestPoint(div_points[-2]))
            pnt2 = cell.GetPoint(locator.FindClosestPoint(div_points[-1]))
            dist1 = math.sqrt(distance(pnt1, div_points[-2]))
            dist2 = math.sqrt(distance(pnt2, div_points[-1]))
            k = -2 if dist1 < dist2 else -1
            O = div_points[k]
            m = rotation_matrix[k+3]

        else:
            m = I
            O = np.array([0, 0, 0])
        
        getData = cell.GetPointData().GetArray(radiusArrayName).GetTuple1
        for j in range(cell.GetNumberOfPoints()):
            point = np.asarray(cell.GetPoints().GetPoint(j))
            tmp = np.dot(np.dot(np.dot(point - O, R), m), R_inv) + O
            centerlinePoints.InsertNextPoint(tmp)
            radiusArray.SetTuple1(count, getData(j))
            centerlineCellArray.InsertCellPoint(count)
            count += 1

    centerline.SetPoints(centerlinePoints)
    centerline.SetLines(centerlineCellArray)
    centerline.GetPointData().AddArray(radiusArray)

    return centerline


def rotationMatrix(data, angle, leave1, leave2):
    """
    Compute the rotation matrices for one or both
    daughter brances of the vessel.

    Args:
        data  (dict): Contains information about landmarking points.
        angle (float): Angle which brances are rotated.
        leave1 (bool): Leaves first daughter branch if True.
        leave2 (bool): Leaves second daughter branch if True.

    Returns:
        R (ndarray): Matrix containing unit vectors in the rotated coordinate system.
    Returns:
        m (dict): Cointains rotation matrices for each daughter branch.
    """

    # Create basis vectors defining bifurcation plane
    d = (np.asarray(data[0]["div_point"]) + \
        np.asarray(data[1]["div_point"]) + \
        np.asarray(data["bif"]["div_point"])) / 3.
    vec = np.eye(3)
    for i in range(2):
        e = np.asarray(data[i]["end_point"])
        tmp = e - d
        len = math.sqrt(np.dot(tmp, tmp))
        vec[:,i] = tmp / len
    
    # Expand basis to 3D
    R = gram_schmidt(vec)
    
    # Set up rotation matrices
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    m1 = np.asarray([[cos_a, -sin_a, 0],
                     [sin_a,  cos_a, 0],
                     [    0,      0, 1]])
    m2 = np.asarray([[ cos_a, sin_a, 0],
                     [-sin_a, cos_a, 0],
                     [     0,     0, 1]])

    m = {1: m1, 2: m2}
    tmp1 = data[0]["div_point"] - d
    tmp2 = data[1]["div_point"] - d

    I = np.eye(3)

    if np.dot(tmp1, R)[0] > np.dot(tmp2, R)[0]:
        m = {1: m2, 2: m1}
    
    # Leave one of the branches untouched
    if leave1:
        k = 1 if data[0]["r_end"] > data[1]["r_end"] else 2
        m[k] = I
    if leave2:
        k = 1 if data[0]["r_end"] < data[1]["r_end"] else 2
        m[k] = I

    return R, m

def get_startpoint(centerline):
    """
    Finds start point of a given centerline.

    Args:
        centerline (vtkPolyData): Centerline data.

    Returns:
        start_point (vtkPoint): Start point of centerline.
    """
    line = extract_single_line(centerline, 0)
    start_point = line.GetPoints().GetPoint(0)
    return start_point

def merge_cl(centerline, end_point, div_point):
    """
    Merge overlapping centerliens.

    Args:
        centerline (vtkPolyData): Centerline data consisting of multiple lines.
        end_point (ndarray): Point where bifurcation ends.
        div_point (ndarray): Point where centerlines diverge. 

    Returns: 
        merge (vtkPolyData): Merged centerline. 
    """
    merge = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    N_lines = centerline.GetNumberOfLines()

    arrays = []
    N_, names = get_number_of_arrays(centerline)
    for i in range(N_):
        tmp = centerline.GetPointData().GetArray(names[i])
        tmp_comp = tmp.GetNumberOfComponents()
        array = get_vtk_array(names[i], tmp_comp, centerline.GetNumberOfPoints())
        arrays.append(array)

    # Find lines to merge
    lines = [extract_single_line(centerline, i) for i in range(N_lines)]
    locators = [get_locator(lines[i]) for i in range(N_lines)]
    div_ID = [locators[i].FindClosestPoint(div_point[0]) for i in range(N_lines)]
    end_ID = [locators[i].FindClosestPoint(end_point[0]) for i in range(N_lines)]
    dist = [np.sum(lines[i].GetPoint(end_ID[i]) - end_point[0]) for i in range(N_lines)]

    # Find the direction of each line
    map_other = {0: 1, 1: 0}
    ID0 = locators[0].FindClosestPoint(end_point[1])
    ID1 = locators[1].FindClosestPoint(end_point[1])
    dist0 = math.sqrt(np.sum((np.asarray(lines[0].GetPoint(ID0)) - end_point[1])**2))
    dist1 = math.sqrt(np.sum((np.asarray(lines[1].GetPoint(ID1)) - end_point[1])**2))
    end1 = 0 if dist0 < dist1 else 1
    end2 = int(not end1)
    for i in range(2, N_lines):
        ID1 = locators[i].FindClosestPoint(end_point[1])
        ID2 = locators[i].FindClosestPoint(end_point[2])
        dist1 = math.sqrt(np.sum((np.asarray(lines[i].GetPoint(ID1)) - end_point[1])**2))
        dist2 = math.sqrt(np.sum((np.asarray(lines[i].GetPoint(ID2)) - end_point[2])**2))
        map_other[i] = end1 if dist1 > dist2 else end2

    counter = 0
    for i in range(centerline.GetNumberOfLines()):
        line = lines[i]

        # Check if it should be merged
        loc = get_locator(line)
        clipp_id = loc.FindClosestPoint(end_point[0])
        div_id = loc.FindClosestPoint(div_point[0])
        clipp_dist = distance(line.GetPoint(clipp_id), end_point[0])
        div_dist = distance(line.GetPoint(div_id), div_point[0])
        tol = get_tolerance(line)*3
        merge_bool = True
        if clipp_dist > tol or div_dist > tol:
            merge_bool = False

        # Get the other line
        other = lines[map_other[i]]
        N = line.GetNumberOfPoints()
        cellArray.InsertNextCell(N)

        for j in range(N):
            # Add point
            if div_ID[i] < j < end_ID[i] and merge_bool: # and i in change:
                new = (np.asarray(other.GetPoint(j)) +
                    np.asarray(line.GetPoint(j))) / 2.
                points.InsertNextPoint(new)
            else:
                points.InsertNextPoint(line.GetPoint(j))

            cellArray.InsertCellPoint(counter)

            # Add array
            for k in range(N_):
                num = arrays[k].GetNumberOfComponents()
                if num == 1:
                    tmp = line.GetPointData().GetArray(names[k]).GetTuple1(j)
                    arrays[k].SetTuple1(counter, tmp)
                elif num == 3:
                    tmp = line.GetPointData().GetArray(names[k]).GetTuple3(j)
                    arrays[k].SetTuple3(counter, tmp[0], tmp[1], tmp[2])
                else:
                    print("Add more options")
                    sys.exit(0)

            counter += 1

    # Insert points, lines and arrays
    merge.SetPoints(points)
    merge.SetLines(cellArray)
    for i in range(N_):
        merge.GetPointData().AddArray(arrays[i])

    return merge

def sort_outlets(outlets, outlet1, outlet2, dirpath):
    """
    Sort all outlets of the geometry given the two relevant outlets

    Args:
        outlets (list): List of outlet center points.
        outlet1 (list): Point representing first relevant oultet. 
        outlet2 (list): Point representing second relevant oultet. 
        dirpath (str): Location of info file.

    Returns:
        outlets (list): List of sorted outlet center points.
    Returns: 
        outlet1 (list): Point representing first relevant oultet. 
    Returns:
        outlet2 (list): Point representing second relevant oultet. 
    """
    tmp_outlets = np.array(outlets).reshape(len(outlets)//3, 3)
    outlet1_index = np.argsort(np.sum((tmp_outlets - outlet1)**2, axis=1))[0]
    outlet2_index = np.argsort(np.sum((tmp_outlets - outlet2)**2, axis=1))[0]
    tmp_outlets = tmp_outlets.tolist()
    if max(outlet1_index, outlet2_index) == outlet1_index:
        outlet1 = tmp_outlets.pop(outlet1_index)
        outlet2 = tmp_outlets.pop(outlet2_index)
    else:
        outlet2 = tmp_outlets.pop(outlet2_index)
        outlet1 = tmp_outlets.pop(outlet1_index)
    outlet_rest = (np.array(tmp_outlets).flatten()).tolist()
    outlets = outlet1 + outlet2 + outlet_rest
    data = {}
    for i in range(len(outlets)//3):
        data["outlet"+str(i)] = outlets[3*i:3*(i+1)]
    write_parameters(data, dirpath)

    return outlets, outlet1, outlet2


def main(dirpath, name, smooth, smooth_factor, angle, l1, l2, bif, lower,
         cylinder_factor, aneurysm, anu_num, resampling_step, version):
    """ 
    Objective rotation of daughter branches, by rotating 
    centerlines and Voronoi diagram about the bifuraction center.
    The implementation is an extension of the original method
    presented by Ford et al. (2009), for aneurysm removal, 
    which introduces the possibility to rotate the 
    daughter branches a given angle. 
    Includes the option to rotate only one of the daughter branches. 

    Args:
        dirpath (str): Directory where case folder is located. 
        name (str): Name of directory where case model is located. 
        smooth (bool): Determine if the voronoi diagram should be smoothed.
        smooth_factor (float): Smoothing factor used for voronoi diagram smoothing.
        angle (float): Angle which daughter branches are moved, in radians. 
        l1 (bool): Leaves first branch untouched if True.
        l2 (bool): Leaves second branch untouched if True.
        bif (bool): Interpolates bifurcation is True.
        lower (bool): Interpolates a lowered line through the bifurcation if True.
        cylinder_factor(float): Factor for choosing the smaller cylinder during Voronoi interpolation.
        aneurysm (bool): Determines if aneurysm is present.
        anu_num (int): Number of aneurysms.
        resampling_step (float): Resampling step used to resample centerlines. 
        version (bool): Determines bifurcation interpolation method.
    """
    # Filenames
    folder = path.join(dirpath, name + "%s")

    # Input filenames
    model_path = folder %  ".vtp"

    # Output names
    # Surface
    model_smoothed_path = folder % "_smooth.vtp"

    # Centerliens
    centerline_par_path = folder % "_centerline_par.vtp"
    centerline_aneurysm_path = folder %  "_centerline_aneurysm.vtp"
    centerline_bif_path = folder %  "_centerline_bif.vtp"
    centerline_complete_path = folder %  "_centerline_complete.vtp"
    centerline_clipped_path = folder %  "_centerline_clipped_ang.vtp"
    centerline_clipped_bif_path = folder %  "_centerline_clipped_bif_ang.vtp"
    centerline_bif_clipped_path = folder %  "centerline_clipped_bif_ang.vtp"
    centerline_dau_clipped_path = folder %  "centerline_clipped_dau_ang.vtp"
    centerline_new_path = folder %  "_centerline_interpolated_ang.vtp"
    centerline_new_bif_path = folder %  "_centerline_interpolated_bif_ang.vtp"
    centerline_new_bif_lower_path = folder %  "_centerline_interpolated_bif_lower_ang.vtp"
    centerline_relevant_outlets_path = folder %  "_centerline_relevant_outlets.vtp"
    centerline_rotated_path = folder %  "centerline_rotated_ang.vtp"
    centerline_rotated_bif_path = folder %  "centerline_rotated_bif_ang.vtp"
    centerline_rotated_dau_path = folder %  "centerline_rotated_dau_ang.vtp"

    # Voronoi diagrams
    voronoi_path = folder % "_voronoi.vtp"
    voronoi_smoothed_path = folder % "_voronoi_smoothed.vtp"
    voronoi_clipped_path = folder % "_voronoi_clipped_ang.vtp"
    voronoi_ang_path = folder %  "_voronoi_ang.vtp"
    voronoi_rotated_path = folder %  "voronoi_rotated_ang.vtp"

    # Points
    points_clipp_path = folder %  "_clippingpoints.vtp"
    points_div_path = folder % "_divergingpoints.vtp"

    # Naming based on different options
    s = "_pi%s" % angle if angle == 0 else "_pi%s" % (math.pi / angle)
    s += "" if not l1 else "_l1"
    s += "" if not l2 else "_l2"
    s += "" if not bif else "_bif"
    s += "" if not smooth else "_smooth"
    s += "" if not lower else "_lower"
    s += "" if cylinder_factor == 7.0 else "_cyl%s" % cylinder_factor
    model_new_surface = folder % ("_angle"+s+".vtp")

    # Read and check model
    if not path.exists(model_path):
        RuntimeError("The given directory: %s did not contain the file: model.vtp" % dirpath)

    # Get aneurysm type
    parameters = get_parameters(dirpath)
    if "aneurysm_type" in parameters.keys():
        aneurysm_type = parameters["aneurysm_type"]
        print("Aneurysm type read from info.txt file: %s" % aneurysm_type)

    # Clean surface
    surface = read_polydata(model_path)
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

   # Check connectivity and only choose the surface with the largest area
    if not "check_surface" in parameters.keys():
        connected_surface = getConnectivity(surface, mode="Largest")
        if connected_surface.GetNumberOfPoints() != surface.GetNumberOfPoints():
            WritePolyData(surface, model_path.replace(".vtp", "_test.vtp"))

    # Get a capped and uncapped version of the surface
    if is_surface_capped(surface):
        open_surface = uncapp_surface(surface)
        capped_surface = surface
    else:
        open_surface = surface
        capped_surface = capp_surface(surface)

    # Get aneurysm "end point"
    if aneurysm:
        aneurysm_point = get_aneurysm_dome(capped_surface, dirpath, anu_num)
    else:
        aneurysm_point = []

    # Get inlet and outlets
    outlet1, outlet2 = get_relevant_outlets(capped_surface, dirpath)
    inlet, outlets = get_centers(open_surface, dirpath)

    # Sort outlets
    outlets, outlet1, outlet2 = sort_outlets(outlets, outlet1, outlet2, dirpath)

    # Compute parent artery and aneurysm centerline
    centerline_par = vmtk_compute_centerlines(inlet, outlets,
                                          centerline_par_path,
                                          capped_surface, resampling=resampling_step)
    centerlines_complete = vmtk_compute_centerlines(inlet, outlets + aneurysm_point,
                                               centerline_complete_path,
                                               capped_surface, resampling=resampling_step)

    # Additional centerline for bifurcation
    centerline_relevant_outlets = vmtk_compute_centerlines(inlet, outlet1 + outlet2,
                                                          centerline_relevant_outlets_path,
                                                          capped_surface,
                                                          resampling=resampling_step)
    centerline_bif = vmtk_compute_centerlines(outlet1, outlet2,
                                         centerline_bif_path,
                                         capped_surface, resampling=resampling_step)

    # Create a tolerance for diverging
    tolerance = get_tolerance(centerline_par)

    # Get data from centerlines and rotation matrix
    data = get_data(centerline_relevant_outlets, centerline_bif, tolerance)
    R, m = rotationMatrix(data, angle, l1, l2) 
    write_parameters(data, dirpath)

    # Compute and smooth voornoi diagram (not aneurysm)
    print("Compute voronoi diagram")
    voronoi = make_voronoi_diagram(surface, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        parameters = get_parameters(dirpath)
        number_of_aneurysms = len([a for a in parameters.keys() if "aneurysm_" in a])
        if number_of_aneurysms == 1:
            voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerline_par, smooth_factor)
        else:
            aneu_centerline = extract_single_line(centerline_complete,
                                                centerline_complete.GetNumberOfCells() - 1)
            div_aneu_id = []
            for i in range(centerline_complete.GetNumberOfCells()-1):
                div_aneu_id.append(centerline_div(aneu_centerline,
                                                  extract_single_line(centerline_complete, i)))
            div_aneu_id = max(div_aneu_id)
            aneu_centerline = extract_single_line(aneu_centerline, start=div_aneu_id)
            voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi,
                                                          centerline_par, smooth_factor,
                                                          no_smooth=aneu_centerline)

        voronoi_smoothed = remove_extreme_points(voronoi_smoothed, voronoi)
        write_polydata(voronoi_smoothed, voronoi_smoothed_path)

        surface_smoothed = create_new_surface(voronoi_smoothed)
        write_polydata(surface_smoothed, model_smoothed_path)

    voronoi = voronoi if not smooth else read_polydata(voronoi_smoothed_path)

    # Locate divpoints and endpoints, for bif or lower, rotated or not
    key = "div_point"
    div_points = get_points(data, key, R, m, rotated=False, bif=False)
    div_points_rotated = get_points(data, key, R, m, rotated=True, bif=False)
    div_points_rotated_bif = get_points(data, key, R, m, rotated=True, bif=True)

    key = "end_point"
    end_points = get_points(data, key, R, m, rotated=False, bif=False)
    end_points_rotated = get_points(data, key, R, m, rotated=True, bif=False)
    end_points_bif = get_points(data, key, R, m, rotated=False, bif=True)
    end_points_rotated_bif = get_points(data, key, R, m, rotated=True, bif=True)

    write_points(div_points[0], points_div_path)
    write_points(end_points[0], points_clipp_path)

    # Clip centerlines 
    print("Clipping centerlines and voronoi diagram.")
    patch_cl = CreateParentArteryPatches(centerline_par, end_points[0])
    write_polydata(patch_cl, centerline_clipped_path)

    if lower or bif:
        patch_bif_cl = CreateParentArteryPatches(centerline_bif, end_points_bif[0])
        write_polydata(patch_bif_cl, centerline_clipped_bif_path)

    # Clip the voronoi diagram
    print("Clipping the Voronoi diagram")
    if path.exists(voronoi_clipped_path):
        voronoi_clipped = read_polydata(voronoi_clipped_path)
    else:
        masked_voronoi = MaskVoronoiDiagram(voronoi, patch_cl)
        voronoi_clipped = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
        write_polydata(voronoi_clipped, voronoi_clipped_path)

    # Rotate branches (Centerline and Voronoi diagram) 
    print("Rotate centerlines and voronoi diagram.")
    rotated_cl = rotate_cl(patch_cl, end_points[1], m, R)
    write_polydata(rotated_cl, centerline_rotated_path)
    
    if lower or bif:
        rotated_bif_cl = rotate_cl(patch_bif_cl, end_points_bif[1], m, R)
        write_polydata(rotated_bif_cl, centerline_rotated_bif_path)
    
    rotated_voronoi = rotate_voronoi(voronoi_clipped, patch_cl, end_points[1], m, R)
    write_polydata(rotated_voronoi, voronoi_rotated_path)


    # Interpolate the centerline
    print("Interpolate centerlines and voronoi diagram.")
    interpolated_cl = InterpolatePatchCenterlines(rotated_cl, centerline_par,
                                                  div_points_rotated[0].GetPoint(0),
                                                  None, False)
    write_polydata(interpolated_cl, centerline_new_path.replace(".vtp", "1.vtp"))


    if bif:
        print("Start interpolate bif")
        interpolated_bif = InterpolatePatchCenterlines(rotated_bif_cl, centerline_bif,
                                                        None, "bif", True)
        write_polydata(interpolated_bif, centerline_new_bif_path)

    if lower:
        print("Start interpolate lower")
        center = ((1/9.)*div_points[1][0] + (4/9.)*div_points[1][1] + \
                        (4/9.)*div_points[1][2]).tolist()
        div_points_rotated_bif[0].SetPoint(0, center[0], center[1], center[2])
        interpolated_bif_lower = InterpolatePatchCenterlines(rotated_bif_cl, centerline_bif,
                                                             div_points_rotated_bif[0].GetPoint(0),
                                                             "lower", True)
        write_polydata(interpolated_bif_lower, centerline_new_bif_lower_path)

    print("Start merge")
    interpolated_cl = merge_cl(interpolated_cl, div_points_rotated[1],
                               end_points_rotated[1])
    write_polydata(interpolated_cl, centerline_new_path)

    bif_ = []
    if lower and bif:
        bif_ = [interpolated_bif, interpolated_bif_lower, rotated_bif_cl]
    elif bif:
        bif_ = [interpolated_bif, rotated_bif_cl]
    elif lower:
        bif_ = [interpolated_bif_lower, rotated_bif_cl]

    # Interpolate voronoi diagram
    print("Start interpolate voronoi diagram")
    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl, rotated_cl,
                                                       rotated_voronoi,
                                                       end_points_rotated,
                                                       bif_, lower, cylinder_factor)

    write_polydata(interpolated_voronoi, voronoi_ang_path.replace(".vtp", "_remove.vtp"))
    interpolated_voronoi = remove_distant_points(interpolated_voronoi, interpolated_cl)
    write_polydata(interpolated_voronoi, voronoi_ang_path)

    # Write a new surface from the new voronoi diagram
    print("Create new surface")
    new_surface = create_new_surface(interpolated_voronoi)

    print("Surface saved in: {}".format(model_new_surface.split("/")[-1]))
    write_polydata(new_surface, model_new_surface)



if  __name__ == "__main__":
    smooth, angle, smooth_factor, l1, l2, bif, basedir, case, lower, \
    cylinder_factor, version, aneurysm, anu_num, resampling_step = read_command_line()
    folders = listdir(basedir) if case is None else [case]
    name = "surface/model"
    for folder in folders:
        if folder[:2] in ["P0", "C0"]:
            print("Working on case", folder)
            main(path.join(basedir, folder), name, smooth, smooth_factor, angle, l1,
                  l2, bif, lower, cylinder_factor, aneurysm, anu_num, resampling_step, version)

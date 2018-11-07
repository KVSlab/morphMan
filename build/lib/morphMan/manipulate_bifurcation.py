##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from argparse_common import *
# Local import
from common import *


def rotate_branches(input_filepath, output_filepath, smooth, smooth_factor, angle,
                    keep_fixed_1, keep_fixed_2, bif, lower, no_smooth, no_smooth_point,
                    poly_ball_size, cylinder_factor, resampling_step,
                    region_of_interest, region_points):
    """
    Objective rotation of daughter branches, by rotating
    centerlines and Voronoi diagram about the bifurcation center.
    The implementation is an extension of the original method
    presented by Ford et al. (2009), for aneurysm removal,
    which introduces the possibility to rotate the
    daughter branches a given angle.
    Includes the option to rotate only one of the daughter branches.

    Args:
        input_filepath (str): Path to input surface.
        output_filepath (str): Path to output surface.
        smooth (bool): Determine if the voronoi diagram should be smoothed.
        smooth_factor (float): Smoothing factor used for voronoi diagram smoothing.
        angle (float): Angle which daughter branches are moved, in radians.
        keep_fixed_1 (bool): Leaves first branch untouched if True.
        keep_fixed_2 (bool): Leaves second branch untouched if True.
        bif (bool): Interpolates bifurcation is True.
        lower (bool): Interpolates a lowered line through the bifurcation if True.
        cylinder_factor(float): Factor for choosing the smaller cylinder during Voronoi interpolation.
        resampling_step (float): Resampling step used to resample centerlines.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' ]
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        poly_ball_size (list): Resolution of polyballs used to create surface.
    """
    # Filenames
    base_path = get_path_names(input_filepath)

    # Output filepaths
    # Centerliens
    centerline_par_path = base_path + "_centerline_par.vtp"
    centerline_bif_path = base_path + "_centerline_bif.vtp"
    centerline_clipped_path = base_path + "_centerline_clipped_ang.vtp"
    centerline_clipped_bif_path = base_path + "_centerline_clipped_bif_ang.vtp"
    centerline_new_path = base_path + "_centerline_interpolated_ang.vtp"
    centerline_new_bif_path = base_path + "_centerline_interpolated_bif_ang.vtp"
    centerline_new_bif_lower_path = base_path + "_centerline_interpolated_bif_lower_ang.vtp"
    centerline_relevant_outlets_path = base_path + "_centerline_relevant_outlets.vtp"
    centerline_rotated_path = base_path + "_centerline_rotated_ang.vtp"
    centerline_rotated_bif_path = base_path + "_centerline_rotated_bif_ang.vtp"
    centerline_bif_clipped_path = base_path + "_centerline_clipped_out.vtp"

    # Voronoi diagrams
    voronoi_clipped_path = base_path + "_voronoi_clipped_ang.vtp"
    voronoi_ang_path = base_path + "_voronoi_ang.vtp"
    voronoi_rotated_path = base_path + "_voronoi_rotated_ang.vtp"

    # Points
    points_clipp_path = base_path + "_clippingpoints.vtp"
    points_div_path = base_path + "_divergingpoints.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Get inlet and outlets
    inlet, outlets = get_centers(surface, base_path)
    if region_of_interest == "manual":
        outlet1, outlet2 = get_relevant_outlets(capped_surface, base_path)
    else:
        outlet1, outlet2 = region_points[:3], region_points[3:]
        surface_locator = get_locator(capped_surface)
        id1 = surface_locator.FindClosestPoint(outlet1)
        id2 = surface_locator.FindClosestPoint(outlet2)
        outlet1 = capped_surface.GetPoint(id1)
        outlet2 = capped_surface.GetPoint(id2)

    print("-- Region of interest is defined by the region points: \nOutlet 1: %s \nOutlet 2: %s" % (outlet1, outlet2))

    # Sort outlets
    outlets, outlet1, outlet2 = sort_outlets(outlets, outlet1, outlet2, base_path)

    # Compute parent artery and aneurysm centerline
    centerline_par, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerline_par_path, capped_surface,
                                                            resampling=resampling_step, base_path=base_path)

    # Additional centerline for bifurcation
    centerline_relevant_outlets, _, _ = compute_centerlines(inlet, outlet1 + outlet2, centerline_relevant_outlets_path,
                                                            capped_surface, resampling=resampling_step, voronoi=voronoi,
                                                            pole_ids=pole_ids, base_path=base_path)

    centerline_bif, _, _ = compute_centerlines(outlet1, outlet2, centerline_bif_path, capped_surface,
                                               resampling=resampling_step, voronoi=voronoi, pole_ids=pole_ids)

    # Create a tolerance for diverging
    tolerance = get_tolerance(centerline_par)

    # Get data from centerlines and rotation matrix
    data = get_data(centerline_relevant_outlets, centerline_bif, tolerance)
    R, m = rotation_matrix(data, angle, keep_fixed_1, keep_fixed_2)
    write_parameters(data, base_path)

    # Compute and smooth voronoi diagram (not aneurysm)
    print("-- Compute voronoi diagram.")
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerline_par, base_path, smooth,
                                          smooth_factor, no_smooth, no_smooth_point,
                                          voronoi, pole_ids)

    # Locate diverging-points and end-points, for bif or lower, rotated or not
    key = "div_point"
    div_points = get_points(data, key, bif=False)

    key = "end_point"
    end_points = get_points(data, key, bif=False)
    end_points_bif = get_points(data, key, bif=True)

    write_points(div_points[0], points_div_path)
    write_points(end_points[0], points_clipp_path)

    # Clip centerlines
    print("-- Clipping centerlines.")
    patch_cl = create_parent_artery_patches(centerline_par, end_points[0])
    write_polydata(patch_cl, centerline_clipped_path)

    # Get the centerline which was clipped away
    clipped_centerline = get_clipped_centerline(centerline_relevant_outlets, data)
    write_polydata(clipped_centerline, centerline_bif_clipped_path)

    if lower or bif:
        patch_bif_cl = create_parent_artery_patches(centerline_bif, end_points_bif[0])
        write_polydata(patch_bif_cl, centerline_clipped_bif_path)

    # Clip the voronoi diagram
    print("-- Clipping the Voronoi diagram")
    voronoi_clipped, _ = split_voronoi_with_centerlines(voronoi, [patch_cl,
                                                                  clipped_centerline])
    write_polydata(voronoi_clipped, voronoi_clipped_path)

    # Rotate branches (Centerline and Voronoi diagram)
    print("-- Rotate centerlines and voronoi diagram.")
    rotated_cl = rotate_cl(patch_cl, end_points[1], m, R)
    write_polydata(rotated_cl, centerline_rotated_path)

    if lower or bif:
        rotated_bif_cl = rotate_cl(patch_bif_cl, end_points_bif[1], m, R)
        write_polydata(rotated_bif_cl, centerline_rotated_bif_path)

    rotated_voronoi = rotate_voronoi(voronoi_clipped, patch_cl, end_points[1], m, R)
    write_polydata(rotated_voronoi, voronoi_rotated_path)

    # Interpolate the centerline
    print("-- Interpolate centerlines.")
    interpolated_cl = interpolate_patch_centerlines(rotated_cl, centerline_par, div_points[0].GetPoint(0),
                                                    None, False)
    write_polydata(interpolated_cl, centerline_new_path.replace(".vtp", "1.vtp"))

    if bif:
        interpolated_bif = interpolate_patch_centerlines(rotated_bif_cl, centerline_bif,
                                                         None, "bif", True)
        write_polydata(interpolated_bif, centerline_new_bif_path)

    if lower:
        center = ((1 / 9.) * div_points[1][0] + (4 / 9.) * div_points[1][1] +
                  (4 / 9.) * div_points[1][2]).tolist()
        div_points[0].SetPoint(0, center[0], center[1], center[2])
        interpolated_bif_lower = interpolate_patch_centerlines(rotated_bif_cl, centerline_bif,
                                                               div_points[0].GetPoint(0),
                                                               "lower", True)
        write_polydata(interpolated_bif_lower, centerline_new_bif_lower_path)

    interpolated_cl = merge_cl(interpolated_cl, div_points[1],
                               end_points[1])
    write_polydata(interpolated_cl, centerline_new_path)

    bif_ = []
    if lower and bif:
        bif_ = [interpolated_bif, interpolated_bif_lower, rotated_bif_cl]
    elif bif:
        bif_ = [interpolated_bif, rotated_bif_cl]
    elif lower:
        bif_ = [interpolated_bif_lower, rotated_bif_cl]

    # Interpolate voronoi diagram
    print("-- Interpolate voronoi diagram.")
    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl, rotated_cl,
                                                       rotated_voronoi,
                                                       end_points,
                                                       bif_, cylinder_factor)
    # Note: This function is slow, and can be commented, but at the cost of robustness.
    interpolated_voronoi = remove_distant_points(interpolated_voronoi, interpolated_cl)
    write_polydata(interpolated_voronoi, voronoi_ang_path)

    # Write a new surface from the new voronoi diagram
    print("-- Create new surface.")
    new_surface = create_new_surface(interpolated_voronoi, poly_ball_size)

    print("-- Preparing surface for output.")
    new_surface = prepare_surface_output(new_surface, surface, interpolated_cl,
                                         output_filepath, test_merge=True, changed=True,
                                         old_centerline=centerline_par)

    print("-- Writing new surface to {}.".format(output_filepath))
    write_polydata(new_surface, output_filepath)


def get_points(data, key, bif=False):
    """
    Finds specific points around the bifurcation, based on the
    key argument. Points can before or after rotation.

    Args:
        data (dict): Contains information about points and IDs of branches and bifurcation.
        key (str): Type of points to extract.
        bif (true): Gets only bifurcation points if True.

    Returns:
        points (vtkPoints): Points as VTK objects.
    Returns:
        div_points_bif (ndarray): Points as numpy objects.
    """
    div_points = np.asarray([data["bif"][key], data[0][key], data[1][key]])

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
        m (dict): Contains rotation matrices for each daughter branch.
    Returns:
        masked_voronoi (vtkPolyData): Rotated voronoi diagram.
    """
    number_of_points = clipped_voronoi.GetNumberOfPoints()
    vtk_distance = vtk.vtkMath.Distance2BetweenPoints
    I = np.eye(3)
    R_inv = np.linalg.inv(R)

    locator = []
    cell_line = []
    not_rotate = [0]
    for i in range(patch_cl.GetNumberOfCells()):
        cell_line.append(extract_single_line(patch_cl, i))
        tmp_locator = get_locator(cell_line[-1])
        locator.append(tmp_locator)

    for i in range(1, patch_cl.GetNumberOfCells()):
        pnt = cell_line[i].GetPoints().GetPoint(0)
        new = cell_line[0].GetPoints().GetPoint(locator[0].FindClosestPoint(pnt))
        dist = math.sqrt(vtk_distance(pnt, new)) < divergingRatioToSpacingTolerance
        if dist:
            not_rotate.append(i)

    def check_rotate(point):
        dist = []
        for i in range(len(locator)):
            tmp_id = locator[i].FindClosestPoint(point)
            tmp = cell_line[i].GetPoints().GetPoint(tmp_id)
            dist.append(math.sqrt(vtk_distance(tmp, point)))

        if dist.index(min(dist)) not in not_rotate:
            pnt = cell_line[dist.index(min(dist))].GetPoints().GetPoint(0)
            if math.sqrt(vtk_distance(pnt, div_points[1])) > \
                    math.sqrt(vtk_distance(pnt, div_points[2])):
                m_ = m[2]
                div = div_points[2]
            else:
                m_ = m[1]
                div = div_points[1]
            return m_, div
        else:
            return I, np.array([0, 0, 0])

    masked_voronoi = vtk.vtkPolyData()
    masked_points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, number_of_points)

    # Iterate through voronoi diagram
    for i in range(number_of_points):
        point = [0.0, 0.0, 0.0]
        clipped_voronoi.GetPoint(i, point)

        point_radius = clipped_voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        M, origo = check_rotate(point)
        tmp = np.dot(np.dot(np.dot(np.asarray(point) - origo, R), M), R_inv) + origo
        masked_points.InsertNextPoint(tmp)
        radius_array.SetTuple1(i, point_radius)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i)

    masked_voronoi.SetPoints(masked_points)
    masked_voronoi.SetVerts(cell_array)
    masked_voronoi.GetPointData().AddArray(radius_array)

    return masked_voronoi


def rotate_cl(patch_cl, div_points, rotation_matrices, R):
    """
    Perform rotation of the centerline representing the
    daughter branches. Rotate along the bifurcation plane
    spanned by two vectors, preserving the angle with
    the rest of the vasculature. Rotation is performed
    using a standard rotational matrix.

    Args:
        patch_cl (vtkPolyData): Clipped centerline representing two daughter branches.
        div_points (ndarray): Contains bifurcation landmarking points.
        rotation_matrices (dict): Contains rotation matrices for each daughter branch.
        R (ndarray): Matrix containing unit vectors in the rotated coordinate system.
    Returns:
        centerline (vtkPolyData): Rotated centerline.
    """
    distance = vtk.vtkMath.Distance2BetweenPoints
    I = np.eye(3)
    R_inv = np.linalg.inv(R)

    number_of_points = patch_cl.GetNumberOfPoints()

    centerline = vtk.vtkPolyData()
    centerline_points = vtk.vtkPoints()
    centerline_cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, number_of_points)

    line0 = extract_single_line(patch_cl, 0)
    locator0 = get_locator(line0)

    # Iterate through points along the centerline
    count = 0
    for i in range(patch_cl.GetNumberOfCells()):
        cell = extract_single_line(patch_cl, i)
        centerline_cell_array.InsertNextCell(cell.GetNumberOfPoints())

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
            origo = div_points[k]
            m = rotation_matrices[k + 3]
        else:
            m = I
            origo = np.array([0, 0, 0])

        radius_array_data = cell.GetPointData().GetArray(radiusArrayName).GetTuple1
        for j in range(cell.GetNumberOfPoints()):
            point = np.asarray(cell.GetPoints().GetPoint(j))
            tmp = np.dot(np.dot(np.dot(point - origo, R), m), R_inv) + origo
            centerline_points.InsertNextPoint(tmp)
            radius_array.SetTuple1(count, radius_array_data(j))
            centerline_cell_array.InsertCellPoint(count)
            count += 1

    centerline.SetPoints(centerline_points)
    centerline.SetLines(centerline_cell_array)
    centerline.GetPointData().AddArray(radius_array)

    return centerline


def rotation_matrix(data, angle, leave1, leave2):
    """
    Compute the rotation matrices for one or both
    daughter branches of the vessel.

    Args:
        data  (dict): Contains information about landmarking points.
        angle (float): Angle which branches are rotated.
        leave1 (bool): Leaves first daughter branch if True.
        leave2 (bool): Leaves second daughter branch if True.

    Returns:
        R (ndarray): Matrix containing unit vectors in the rotated coordinate system.
    Returns:
        m (dict): Contains rotation matrices for each daughter branch.
    """

    # Create basis vectors defining bifurcation plane
    d = (np.asarray(data[0]["div_point"]) +
         np.asarray(data[1]["div_point"]) +
         np.asarray(data["bif"]["div_point"])) / 3.
    vec = np.eye(3)
    for i in range(2):
        e = np.asarray(data[i]["end_point"])
        tmp = e - d
        length = math.sqrt(np.dot(tmp, tmp))
        vec[:, i] = tmp / length

    # Expand basis to 3D
    R = gram_schmidt(vec)

    # Set up rotation matrices
    cos_a = math.cos(angle)
    sin_a = math.sin(angle)
    m1 = np.asarray([[cos_a, -sin_a, 0],
                     [sin_a, cos_a, 0],
                     [0, 0, 1]])
    m2 = np.asarray([[cos_a, sin_a, 0],
                     [-sin_a, cos_a, 0],
                     [0, 0, 1]])

    m = {1: m1, 2: m2}
    tmp1 = data[0]["div_point"] - d
    tmp2 = data[1]["div_point"] - d

    I = np.eye(3)

    if np.dot(tmp1, R)[0] > np.dot(tmp2, R)[0]:
        m = {1: m2, 2: m1}

    # Leave one of the branches untouched
    if leave1:
        k = 1
        m[k] = I
    if leave2:
        k = 2
        m[k] = I

    return R, m


def merge_cl(centerline, end_point, div_point):
    """
    Merge overlapping centerlines.

    Args:
        centerline (vtkPolyData): Centerline data consisting of multiple lines.
        end_point (ndarray): Point where bifurcation ends.
        div_point (ndarray): Point where centerlines diverge.

    Returns:
        merge (vtkPolyData): Merged centerline.
    """
    merge = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
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

    # Find the direction of each line
    map_other = {0: 1, 1: 0}
    ID0 = locators[0].FindClosestPoint(end_point[1])
    ID1 = locators[1].FindClosestPoint(end_point[1])
    dist0 = math.sqrt(np.sum((np.asarray(lines[0].GetPoint(ID0)) - end_point[1]) ** 2))
    dist1 = math.sqrt(np.sum((np.asarray(lines[1].GetPoint(ID1)) - end_point[1]) ** 2))
    end1 = 0 if dist0 < dist1 else 1
    end2 = int(not end1)
    for i in range(2, N_lines):
        ID1 = locators[i].FindClosestPoint(end_point[1])
        ID2 = locators[i].FindClosestPoint(end_point[2])
        dist1 = math.sqrt(np.sum((np.asarray(lines[i].GetPoint(ID1)) - end_point[1]) ** 2))
        dist2 = math.sqrt(np.sum((np.asarray(lines[i].GetPoint(ID2)) - end_point[2]) ** 2))
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
        tol = get_tolerance(line) * 3
        merge_bool = True
        if clipp_dist > tol or div_dist > tol:
            merge_bool = False

        # Get the other line
        other = lines[map_other[i]]
        N = line.GetNumberOfPoints()
        cell_array.InsertNextCell(N)

        for j in range(N):
            # Add point
            if div_ID[i] < j < end_ID[i] and merge_bool:
                new = (np.asarray(other.GetPoint(j)) +
                       np.asarray(line.GetPoint(j))) / 2.
                points.InsertNextPoint(new)
            else:
                points.InsertNextPoint(line.GetPoint(j))

            cell_array.InsertCellPoint(counter)

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
                    print("-- Add more options")
                    sys.exit(0)

            counter += 1

    # Insert points, lines and arrays
    merge.SetPoints(points)
    merge.SetLines(cell_array)
    for i in range(N_):
        merge.GetPointData().AddArray(arrays[i])

    return merge


def read_command_line():
    """
    Read arguments from commandline
    """

    description = "Removes the bifurcation (possibly with an aneurysm), after which the" + \
                  " daughter branches can be rotated."
    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    # Add common arguments
    add_common_arguments(parser)

    # Set region of interest:
    parser.add_argument("-r", "--region-of-interest", type=str, default="manual",
                        choices=["manual", "commandline"],
                        help="The method for defining the region to be changed. There are" +
                             " two options: 'manual' and 'commandline'. In" +
                             " 'manual' the user will be provided with a visualization of the" +
                             " input surface, and asked to provide an end and start point of the" +
                             " region of interest. Note that not all algorithms are robust over" +
                             " bifurcations. If 'commandline' is provided, then '--region-points'" +
                             " is expected to be provided.")
    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" +
                             " argument have to be given. The method expects two points" +
                             " which defines the start and end of the region of interest. If" +
                             " 'method' is set to stenosis, then one point can be provided as well," +
                             " which is assumed to be the center of a new stenosis." +
                             " Example providing the points (1, 5, -1) and (2, -4, 3):" +
                             " --stenosis-points 1 5 -1 2 -4 3")

    # Arguments for rotation
    parser.add_argument('-a', '--angle', type=float, default=10,
                        help="Each daughter branch is rotated an angle 'a' in the" +
                             " bifurcation plane. 'a' is assumed to be in degrees," +
                             " and not radians", metavar="rotation_angle")
    parser.add_argument("--keep-fixed-1", type=str2bool, default=False,
                        help="Leave one branch untouched")
    parser.add_argument("--keep-fixed-2", type=str2bool, default=False,
                        help="Leave one branch untouched")

    # Bifurcation reconstruction arguments
    parser.add_argument("--bif", type=str2bool, default=False,
                        help="interpolate bif as well")
    parser.add_argument("--lower", type=str2bool, default=False,
                        help="Make a fourth line to interpolate along that" +
                             " is lower than the other bif line.")
    parser.add_argument("--cylinder-factor", type=float, default=7.0,
                        help="Factor for choosing the smaller cylinder")

    args = parser.parse_args()
    ang_ = 0 if args.angle == 0 else args.angle * math.pi / 180  # Convert from deg to rad

    return dict(input_filepath=args.ifile, smooth=args.smooth, output_filepath=args.ofile,
                smooth_factor=args.smooth_factor, angle=ang_,
                keep_fixed_1=args.keep_fixed_1, keep_fixed_2=args.keep_fixed_2,
                bif=args.bif, lower=args.lower, cylinder_factor=args.cylinder_factor,
                resampling_step=args.resampling_step, no_smooth=args.no_smooth,
                no_smooth_point=args.no_smooth_point, poly_ball_size=args.poly_ball_size,
                region_of_interest=args.region_of_interest,
                region_points=args.region_points)


if __name__ == "__main__":
    rotate_branches(**read_command_line())

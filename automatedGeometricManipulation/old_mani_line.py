from common import *
from argparse import ArgumentParser
from os import path, listdir
from subprocess import STDOUT, check_output
from IPython import embed

import operator
import sys
import math
import numpy.linalg as la
import matplotlib.pyplot as plt

# Local import
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *

NUM = 20

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('-d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('-case', type=str, default=None, help="Choose case")
    parser.add_argument('-s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram (surface) should be" + \
                        "smoothed before it is manipulated", metavar="smooth")
    parser.add_argument("-a", "--alpha", type=float, default=None,
                        help="Compression factor in vertical direction, ranging from -1.0 to 1.0")
    parser.add_argument("-b", "--beta", type=float, default=None,
                        help="Compression factor in horizontal direction, ranging from -1.0 to 1.0")

    args = parser.parse_args()

    return args.smooth, args.dir_path, args.case, args.alpha, args.beta


def get_points(data, key, bif=False):
    if bif:
        div_points = np.asarray([data[0][key], data[1][key]])
    else:
        div_points = np.asarray([data["bif"][key], data[0][key], data[1][key]])

    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)

    return points, div_points


def get_clipping_points(dirpath):
    particles = path.join(dirpath, "optimal.particles")
    all_points = np.loadtxt(particles)
    if dirpath[-5:] == "P0252":
        clipping_points = all_points[1::2] 
    else:
        clipping_points = all_points[2:] 
    return clipping_points

def reconst_and_move_voro(dirpath, smooth, name, alpha=0.0, beta=0.0):
    case = dirpath[-5:]
    # Input filenames
    model_path = path.join(dirpath, name, "model.vtp")

    # Output names
    model_smoothed_path = path.join(dirpath, name, "model_smoothed.vtp")
    model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
    model_new_surface_tmp = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s_tmp.vtp" % (alpha, beta))
    print dirpath
    new_centerlines_path = path.join("./modelscurv", "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    model_new_surface = path.join("./modelscurv", "%s_alpha_%s_beta_%s.vtp" % (case, alpha,beta))

    # Centerlines
    centerline_complete_path = path.join(dirpath, name, "centerline_complete.vtp")
    centerline_clipped_path = path.join(dirpath, name,  "centerline_clipped.vtp")
    centerline_clipped_part_path = path.join(dirpath, name,  "centerline_clipped_part.vtp")
    centerline_clipped_part1_path = path.join(dirpath, name,  "centerline_clipped_end_part.vtp")
    centerline_new_path = path.join(dirpath, name, "centerline_interpolated.vtp")

    # Voronoi diagrams
    voronoi_path = path.join(dirpath, name, "model_voronoi.vtp")
    voronoi_smoothed_path = path.join(dirpath, name,  "model_voronoi_smoothed.vtp")
    voronoi_clipped_path = path.join(dirpath, name, "model_voronoi_clipped.vtp")
    voronoi_clipped_part_path = path.join(dirpath, name, "model_voronoi_clipped_part.vtp")
    
    # Extract Clipping points
    clipping_points = get_clipping_points(dirpath)

    # Read and check model
    if not path.exists(model_path):
        RuntimeError("The given directory: %s did not contain the file: model.vtp" % dirpath)

    # Clean surface
    surface = ReadPolyData(model_path)
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

    #Get a capped and uncapped version of the surface
    open_surface = surface
    capped_surface = capp_surface(surface)

    # Get inlet and outlets
    inlet, outlets = get_centers(open_surface, dirpath)
    
    # Compute all centerlines
    centerlines_complete = compute_centerlines(inlet, outlets,
                                               centerline_complete_path,
                                               capped_surface, resampling=0.1)

    print("Compute voronoi diagram")
    voronoi = makeVoronoiDiagram(surface, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerlines_complete, 0.25)
        WritePolyData(voronoi_smoothed, voronoi_smoothed_path)
        surface_smoothed = create_new_surface(voronoi_smoothed)
        WritePolyData(surface_smoothed, model_smoothed_path)

    voronoi = voronoi if not smooth else ReadPolyData(voronoi_smoothed_path)

    # Set clipping points
    div_points = np.asarray(clipping_points)
    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)
    clip_points = points

    # Sort centerlines 
    lines = []
    n = centerlines_complete.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines_complete, i))

    longest = [lines[0]]
    longest_length = get_curvilinear_coordinate(longest[0])

    for i in range(1,n):
        tmplong = get_curvilinear_coordinate(lines[i])
        if len(tmplong) > len(longest_length):
            longest_length = tmplong
            longest.insert(0, lines[i])
        else: 
            longest.append(lines[i])

    centerlines_in_order = merge_data(longest)


    # Check if case includs the opthalmic artery
    eye = False
    eye, clip_ID, centerlines_in_order, eyeline = find_ophthalmic_artery(centerlines_in_order, clipping_points)

    if eye == True:
        print("Clipping opthamlic artery")
        patch_eye = clip_eyeline(eyeline, clipping_points[0], clip_ID)

    # Clip centerline
    print("Clipping centerlines.")
    patch_cl = CreateParentArteryPatches(centerlines_in_order, clip_points)
    clipped_curve = get_clipped_curve(centerlines_in_order, clipping_points) 
    WritePolyData(patch_cl, centerline_clipped_path)
    
    if eye == True:
        eyeline_end = ExtractSingleLine(patch_eye,1)
        clipped_curve = merge_data([clipped_curve, eyeline_end])

    WritePolyData(clipped_curve, centerline_clipped_part_path )


    # Clipp the voronoi diagram
    # And keep the clipped part
    print("Computing clipped voronoi diagrams")
    if not path.exists(voronoi_clipped_path):
        masked_voronoi_clip  = MaskVoronoiDiagram(voronoi, clipped_curve)
        voronoi_clipped_part = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi_clip)
        WritePolyData(voronoi_clipped_part, voronoi_clipped_part_path)
    else:
        voronoi_clipped_part = ReadPolyData(voronoi_clipped_part_path)

    if not path.exists(voronoi_clipped_part_path):
        masked_voronoi  = MaskVoronoiDiagram(voronoi, patch_cl)
        voronoi_clipped = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
        WritePolyData(voronoi_clipped, voronoi_clipped_path)
    else:
        voronoi_clipped = ReadPolyData(voronoi_clipped_path)
    
    # Extract translation vectors
    print("Computing translation directions.")
    line = ExtractSingleLine(centerlines_in_order, 0)
    growDir = "horizont"
    middle_points, middleIds = get_spline_points(line, beta, growDir, case, clip_points,voronoi)

    # Move anterior bend horizontally 
    locator = get_locator(line)
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)
    dx_p1 = middle_points[0] - p1
    dx_p2 = middle_points[-1] - p2
    middle_points = middle_points[1:-1]
    newpoints = []

    if beta != 0.0:
        # Iterate over points P from Voronoi diagram,
        # and move them 
        print("Adjusting voronoi diagram")
        voronoi_clipped = move_clipped_voro(dx_p1, dx_p2, voronoi_clipped, patch_cl, ID1,ID2, clip_ID, clip=False)   
        voronoi_clipped_part = move_clipped_voro(dx_p1, dx_p2, voronoi_clipped_part, clipped_curve, ID1,ID2, clip_ID, clip=True, eye=eye)   
        newVoronoi = merge_data([voronoi_clipped,voronoi_clipped_part])
        new_surface = create_new_surface(newVoronoi)

        # Move line manually for postprocessing
        # Update inlet and outlet positions
        print("Adjusting line manually")
        line1 = ExtractSingleLine(patch_cl, 0)
        lines = []

        for i in range(n-1):
            lines.append(ExtractSingleLine(patch_cl, i+1))
            
        patch_cl_start = move_line_horo(line1, ID1, ID2, dx_p1, clip=False, side="left")
        patch_cl_ends = []
        for line_n in lines:
            patch_cl_end = move_line_horo(line_n, ID1, ID2, dx_p1, clip=False, side="right")
            newpoints.append(patch_cl_end.GetPoint(patch_cl_end.GetNumberOfCells() - 1))
            patch_cl_ends.append(patch_cl_end)


        clipped_part_new = move_line_horo(clipped_curve, ID1, ID2, dx_p1, clip=True, eye=eye)

        if eye:
            newpoints.append(clipped_part_new.GetPoint(clipped_part_new.GetNumberOfCells() - 1))
        
        newpoints.append(patch_cl_start.GetPoint(0))


    else:
        newVoronoi = voronoi
        new_surface = surface

    WritePolyData(new_surface, model_new_surface_tmp)

    if alpha != 0.0:
        if beta == 0.0:
            print("No horizontal movement. Initiating vertical movement")
            new_centerline = centerlines_complete
            WritePolyData(new_centerline, new_centerlines_path) 
        else:
            print "Creating new_centerline_complete.vtp of horizontally moved model"
            new_centerline = makeCenterline(model_new_surface_tmp, new_centerlines_path, smooth=False, resampling=False, newpoints=newpoints, recompute=True)
            WritePolyData(new_centerline, new_centerlines_path) 

        # Vertical movement
        print("Moving geometry vertically")
        translate_vertically(new_surface, new_centerline, clip_points, alpha , dirpath, newVoronoi, voronoi_clipped, voronoi_clipped_part, eye, newpoints, case)
    
    else:
        # Write a new surface from the new voronoi diagram
        WritePolyData(new_surface, model_new_surface)
    return new_surface

def connect_line(line):
    
    # Create edges between new_centerline points
    pts = vtk.vtkPoints()
    for i in range((line.GetNumberOfCells())):
        pts.InsertNextPoint(line.GetPoint(i))

    lines = vtk.vtkCellArray()
    for i in range((line.GetNumberOfCells()) - 1):
        newline = vtk.vtkLine()
        newline.GetPointIds().SetId(0, i)
        newline.GetPointIds().SetId(1, i + 1)
        lines.InsertNextCell(newline)

    line = vtk.vtkPolyData()
    line.SetPoints(pts)
    line.SetLines(lines)
    return line

def move_line_horo(patch_cl, ID1, ID2, dx_p1, clip=False, eye=False, side=None):

    centerline_loc = get_locator(patch_cl)
    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if clip == True:
        if eye == True:
            l1 = ExtractSingleLine(patch_cl, 0)
            l2 = ExtractSingleLine(patch_cl, 1)
            l3 = ExtractSingleLine(patch_cl, 2)
            test_cl = merge_data([l1,l2,l3])

            ID1 = 0
            ID2 = len(get_curvilinear_coordinate(test_cl))
            idmid = int( (ID1 + ID2)/2.)
        else:
            ID1 = 0
            ID2 = len(get_curvilinear_coordinate(patch_cl))
            idmid = int( (ID1 + ID2)/2.)
        
        for p in range(patch_cl.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(patch_cl.GetPoint(p))

            if cl_id < idmid:
                dist = dx_p1 * (idmid**2-cl_id**2) / (idmid**2-ID1**2)
            else:
                if cl_id <= (ID2 - 1):
                    dist = -dx_p1 * (cl_id-idmid)**(0.5) / (ID2-idmid)**(0.5)
                else: 
                    locator = get_locator(test_cl)
                    pp = patch_cl.GetPoint(cl_id)
                    id_main = locator.FindClosestPoint(pp)
                    dist = -dx_p1 * (id_main-idmid)**(0.5) / (id_main-idmid)**(0.5)

            patch_point = np.asarray(patch_cl.GetPoint(p))

            if la.norm(patch_point) > 0.1:
                points.InsertNextPoint(patch_point + dist)
                verts.InsertNextCell(1)
                verts.InsertCellPoint(p)

    else:
        for p in range(patch_cl.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(patch_cl.GetPoint(p))

            if side == "right":
                dist = -dx_p1
            elif side == "left":
                dist = dx_p1


            points.InsertNextPoint(np.asarray(patch_cl.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(patch_cl.GetPointData().GetArray(radiusArrayName))

    return newline

def spline_line(line, get_curv=False, isline=False):
    
    if isline == False:
        # Create edges between new_centerline points
        pts = vtk.vtkPoints()
        for i in range((line.GetNumberOfCells())):
            pts.InsertNextPoint(line.GetPoint(i))

        lines = vtk.vtkCellArray()
        for i in range((line.GetNumberOfCells()) - 1):
            newline = vtk.vtkLine()
            newline.GetPointIds().SetId(0, i)
            newline.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(newline)

        line = vtk.vtkPolyData()
        line.SetPoints(pts)
        line.SetLines(lines)

    # Collect data from centerline
    data = np.zeros((line.GetNumberOfPoints(), 3))

    curv_coor = get_curvilinear_coordinate(line)
    for i in range(data.shape[0]):
        data[i,:] = line.GetPoint(i)

    nknots = 50
    t = np.linspace(curv_coor[0], curv_coor[-1], nknots+2)[1:-1]
    fx = splrep(curv_coor, data[:,0], k=4, t=t)
    fy = splrep(curv_coor, data[:,1], k=4, t=t)
    fz = splrep(curv_coor, data[:,2], k=4, t=t)

    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    data = np.zeros((len(curv_coor), 3))
    data[:,0] = fx_
    data[:,1] = fy_
    data[:,2] = fz_

    header = ["X", "Y", "Z"]
    line = data_to_vtkPolyData(data, header)
   
    # Let vmtk compute curve attributes
    line = CenterlineAttribiutes(line, smooth=False)

    if get_curv == True:
        # Analytical curvature
        dlsfx = splev(curv_coor, fx, der=1)
        dlsfy = splev(curv_coor, fy, der=1)
        dlsfz = splev(curv_coor, fz, der=1)

        ddlsfx = splev(curv_coor, fx, der=2)
        ddlsfy = splev(curv_coor, fy, der=2)
        ddlsfz = splev(curv_coor, fz, der=2)

        C1xC2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
        C1xC2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
        C1xC2_3 = ddlsfy * dlsfx - ddlsfx * dlsfy

        curvature = np.sqrt(C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2) / \
                            (dlsfx**2 + dlsfy**2 + dlsfz**2)**1.5

        return line, curvature
    else:
        return line


def translate_vertically(surface, centerline, clip_points, alpha, dirpath, voronoi, voronoi_clipped, voronoi_clipped_part, eye, oldpoints, case):
    # Filenames
    #new_centerlines_path = path.join(dirpath, name, "new_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta))
    #model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))

    new_centerlines_path = path.join("./modelscurv", "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    model_new_surface = path.join("./modelscurv", "%s_alpha_%s_beta_%s.vtp" % (case, alpha,beta))

    # Vertical movement
    locator = get_locator(centerline)
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)
    p1      = centerline.GetPoint(ID1)
    p2      = centerline.GetPoint(ID2)

    # Create clipping points as VTK points
    clipping_points = [p1,p2]
    div_points = np.asarray(clipping_points)
    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)
    clip_points = points

    lines = []
    n = centerline.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerline, i))

    longest = [lines[0]]
    lenlong = get_curvilinear_coordinate(longest[0])

    # Swap longest with first element
    for i in range(1,n):
        tmplong = get_curvilinear_coordinate(lines[i])
        if len(tmplong) > len(lenlong):
            lenlong = tmplong
            longest.insert(0, lines[i])
        else: 
            longest.append(lines[i])

    centerlines_in_order = merge_data(longest)

    # Special cases including the ophthalmic artery
    clip_ID = None
    if eye == True:
        eye, clip_ID, centerlines_in_order, eyeline = find_ophthalmic_artery(centerlines_in_order, clipping_points)

    if eye == True:
        print("Clipping opthamlic artery")
        patch_eye = clip_eyeline(eyeline, clipping_points[0], clip_ID)

    
    print("Clipping centerlines.")
    patch_cl = CreateParentArteryPatches(centerlines_in_order, clip_points)
    
    # Get clipped curve 
    print "Getting clipped curve"
    clipped_curve = get_clipped_curve(centerlines_in_order, clipping_points) 

    if eye == True:
        eyeline_end = ExtractSingleLine(patch_eye,1)
        clipped_curve = merge_data([clipped_curve, eyeline_end])

    # Find ID of middle pooint:
    print("Finding points to spline through.")
    line = ExtractSingleLine(centerlines_in_order, 0)

    loc = get_locator(line)
    ID1 = loc.FindClosestPoint(p1)

    growDir = "vertical"
    middle_points, middleIds, dx = get_spline_points(line, alpha, growDir, case, clip_points,voronoi)
    # Iterate over points P from Voronoi diagram,
    # and move them 
    print("Adjust voronoi diagram")
    moved_clip =  move_voro(voronoi_clipped_part,clipped_curve, ID1,  clip_ID, dx,  eye)
    newVoronoi = merge_data([voronoi_clipped, moved_clip]) 


    # Move centerline manually for postprocessing  
    print("Adjusting line manually")
    if eye:
        newpoints = oldpoints[:-2] + [oldpoints[-1]]
    else:
        newpoints = oldpoints

    clipped_part_new = move_dots_vertical(clipped_curve, dx, eye=eye)

    if eye:
        newpoints.insert(0, clipped_part_new.GetPoint(clipped_part_new.GetNumberOfCells() - 1))

    


    # Write a new surface from the new voronoi diagram
    print("Writing new surface")
    new_surface = create_new_surface(newVoronoi)
    WritePolyData(new_surface, model_new_surface)

    model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
    WritePolyData(new_surface, model_new_surface)

    print "Creating new_centerline_complete.vtp of vertically moved model"
    model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
    new_centerline = makeCenterline(model_new_surface, new_centerlines_path, smooth=False, resampling=False, newpoints=newpoints, recompute=True)
    WritePolyData(new_centerline, new_centerlines_path) 

def move_dots_vertical(line, dx,eye=False):

    centerline_loc = get_locator(line)

    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if eye == True:
        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(line))#test_cl))
        IDmid = int( (ID1 + ID2)/2.)

        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))

            if cl_id <= (ID2-1):
                dist = 4 * dx * (cl_id - ID1)*(ID2 - cl_id) / ( ID2 - ID1)**2
            else:
                dist = dx

            points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    else:
        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(line))
    
        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))

            dist = 4 * dx * (cl_id - ID1)*(ID2 - cl_id) / ( ID2 - ID1)**2

            points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(line.GetPointData().GetArray(radiusArrayName))
    return newline

def find_ophthalmic_artery(centerlines, clip_pts):

    # Extract lines:
    lines = []
    n = centerlines.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines, i))
        WritePolyData(lines[i], "line%i.vtp"%i)

    longest = lines[0]
    tol = 0.40
    
    # Find start and stop IDs along clipped curve
    locator_longest = get_locator(longest)
    p1      = clip_pts[0]
    p2      = clip_pts[1]
    ID1     = locator_longest.FindClosestPoint(p1)
    ID2     = locator_longest.FindClosestPoint(p2)
    
    eye = False
    index = 1
    for line in lines[1:]:
        locator_checkline = get_locator(line)
        len_check = len(get_curvilinear_coordinate(line)) 
        if len_check < ID2:
            IDStop = len_check - 1
        else:
            IDStop = ID2

        for i in np.arange(ID1, IDStop):
            p_eye   = np.asarray(line.GetPoint(i))
            p_cl    = np.asarray(longest.GetPoint(i))
            dist = la.norm(p_eye - p_cl)
            if  dist > tol:
                clip_ID = i
                centerlines = merge_data(lines)
                eye = True
                eyeline = lines[index]
                del lines[index]
                centerlines = merge_data(lines)
                break
        if eye:
            break
        index += 1
    
    if eye:
        return eye, clip_ID, centerlines, eyeline
    else:
        return eye, None, centerlines, None


def get_spline_points(line, param, growDir, case, clip_points, voronoi):
    """ 
    Pick n uniformly selected points along the 
    centerline from point P1 to P2
    And move them
    """
    locator = get_locator(line)
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)
    ID_mid = int((ID1 + ID2) / 2.)
    P = [p1,p2]
    
    # Uniform selection of points along siphon
    points = []
    n = 7
    ids = np.zeros(n)
    dx = 1 / (n + 1.)
    for i in range(1,n+1):
        ID = int(ID1 + (ID2 - ID1) * i * dx)
        ids[i-1] = ID
        p = line.GetPoints().GetPoint(ID)
        points.append(np.array([p[0], p[1], p[2]]))

    for i in range(len(P)):
        P[i] = np.array([P[i][0], P[i][1], P[i][2]])

    xx,yy,zz, n = bestPlane(points,P)

    if growDir == "vertical":
        dz, dx = movePerp(xx,yy,zz,n,P, points, param)
        return dz, ids, dx
    else: 
        dz,zp,zm = movePara(xx,yy,zz,n,P, points, param, case)
        return dz, ids

def bestPlane(Z, P):
    # Defined matrices
    A = Z
    C = P
    b = np.ones(len(A))
    d = np.ones(len(C))

    # Create complete matrix
    ATA = np.transpose(A).dot(A)
    M0 = np.c_[ATA,np.transpose(C)]
    M1 = np.c_[C, np.zeros( (len(C), len(C)))]
    M = np.r_[M0, M1]
    Y = np.r_[np.transpose(A).dot(b),d]
    
    # Solve system
    x = la.solve(M,Y)
    a = x[0]
    b = x[1]
    c = x[2]
    n = np.array([a,b,c])
    n = n / la.norm(n)

    # Define plane
    xmin = min(Z, key=operator.itemgetter(1))[0] - 4 
    xmax = max(Z, key=operator.itemgetter(1))[0] + 4
    ymin = min(Z, key=operator.itemgetter(1))[1] - 4 
    ymax = max(Z, key=operator.itemgetter(1))[1] + 4
    xx,yy = np.meshgrid(np.linspace(xmin,xmax,15),np.linspace(ymin,ymax,15))
    zz = (1 - a*xx - b*yy ) / float(c)
    return xx,yy,zz,n

def movePerp(xx,yy,zz,n, P, Z, alpha):
    p1 = P[0]
    p2 = P[1]

    # Find midpoint and point furthest away
    dist = []
    for z in Z:
        d = la.norm( np.cross((z - p1),(z - p2))) / la.norm(p2 - p1)
        dist.append(d)

    D_id = dist.index(max(dist))
    D = max(dist)
    z_m = Z[D_id]

    # Vector from line to Z_max and projection onto plane
    v = (z_m - p2) - (z_m - p2).dot(p2-p1)*(p2-p1) / la.norm(p2-p1)**2
    PV = v - v.dot(n)*n

    # Find distances 
    P1 =  (z_m - p1).dot(p2-p1)*(p2-p1) / la.norm(p2-p1)**2
    P1 = p1 + P1
    V = P1 + v
    PV1 = P1 + PV

    # Move points
    dZ = []
    for i in range(len(Z)):
        dz = np.array(PV) * dist[i] / D * alpha
        dZ.append( Z[i] + dz)

    return dZ, (PV1 - P1)*alpha

def movePara(xx,yy,zz,n, P,Z, beta, case ):
    p1 = P[0]
    p2 = P[1]
    
    # Find normal from q
    q = [(p1 + p2)/2.]
    qp2 = np.array(p2 - q)[0]
    qp1 = np.array(p1 - q)[0]
    s = q[0] - np.cross(qp2, n)*3

    # Split points based on orientation
    # to q normal
    Z_p = []
    Z_p_dist = []
    Z_m = []
    Z_m_dist = []
    for z in Z:
        d = la.norm( np.cross((z - s),(z - q[0]))) / la.norm(q[0] - s)
        c = np.cross(s-q, z-q)
        if c[0][0] >= 0:
            Z_p.append(z)
            Z_p_dist.append(d)
        else:
            Z_m.append(z)
            Z_m_dist.append(d)
    
    # Move points
    D = la.norm(p1-p2) / 2.
    dZ = []
    dZ.append(p1 + qp1*beta)
    for i in range(len(Z_p)):
        dz = qp1 * Z_p_dist[i] / D * beta
        dZ.append(Z_p[i] + dz)

    for i in range(len(Z_m)):
        dz = qp2 * Z_m_dist[i] / D * beta
        dZ.append(Z_m[i] + dz)
    dZ.append(p2 + qp2*beta)
    
    # Check if moved in right direction
    d_0 = la.norm( np.cross(Z_p[0] - q, Z_p[0] - s) ) / la.norm(s-q)
    d_1 = la.norm( np.cross(dZ[1] - q, dZ[1] - s) ) / la.norm(s-q)
    if d_1 < d_0:
        # Split points based on orientation
        # to q normal
        Z_p = []
        Z_p_dist = []
        Z_m = []
        Z_m_dist = []
        for z in Z:
            d = la.norm( np.cross((z - s),(z - q[0]))) / la.norm(q[0] - s)
            c = -np.cross(s-q, z-q)
            if c[0][0] >= 0:
                Z_p.append(z)
                Z_p_dist.append(d)
            else:
                Z_m.append(z)
                Z_m_dist.append(d)

        # Move points
        D = la.norm(p1-p2) / 2.
        dZ = []
        dZ.append(p1 + qp1*beta)
        for i in range(len(Z_p)):
            dz = qp1 * Z_p_dist[i] / D * beta
            dZ.append(Z_p[i] + dz)

        for i in range(len(Z_m)):
            dz = qp2 * Z_m_dist[i] / D * beta
            dZ.append(Z_m[i] + dz)
        dZ.append(p2 + qp2*beta)


    zpid = Z_p_dist.index(min(Z_p_dist))
    zp_min = Z_p[zpid]
    zmid = Z_m_dist.index(min(Z_m_dist))
    zm_min = Z_m[zmid]
    return dZ, zp_min, zm_min

def get_clipped_curve(centerlines_in_order, clipping_points, end=False): 
    # Extract clipped segment
    start_p = vtk.vtkPoints()
    end_p   = vtk.vtkPoints()
    line    = ExtractSingleLine(centerlines_in_order, 0)
    start   = line.GetPoint(0)
    endId   = len(get_curvilinear_coordinate(line)) -2
    end     = line.GetPoint(endId)
    for p1, p2 in zip([start, clipping_points[0]], [clipping_points[1], end]):
        start_p.InsertNextPoint(p1)
        end_p.InsertNextPoint(p2)

    clip_end = CreateParentArteryPatches(line, end_p,True)

    if end == True:
        cut_end = vtk.vtkPoints()
        cut_end.InsertNextPoint(start)
        cut_end.InsertNextPoint(clipping_points[1])
        return CreateParentArteryPatches(centerlines_in_order, cut_end, True)
    else: 
        clipped_curve = CreateParentArteryPatches(clip_end, start_p,True)

        return clipped_curve

def clip_eyeline(eyeline, clip_start_point, clip_end_ID):

    points = [clip_start_point, eyeline.GetPoint(clip_end_ID)]
    eye_points = vtk.vtkPoints()
    for p in points:
        eye_points.InsertNextPoint(p)

    patch_eye = CreateParentArteryPatches(eyeline, eye_points) 
    return patch_eye

def move_clipped_voro(dx_p1, dx_p2, voronoi_clipped, patch_cl, ID1,ID2, clip_ID, clip=False, eye=False):   

    centerline_loc = get_locator(patch_cl)
    newDataSet = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()
    idmid_0 = int( (ID1 + ID2) / 2.)

    if clip == True:
        ID1_0 = ID1
        ID1 = 0
        if eye == True:
            l1 = ExtractSingleLine(patch_cl, 0)
            l2 = ExtractSingleLine(patch_cl, 1)
            l3 = ExtractSingleLine(patch_cl, 2)
            test_cl = merge_data([l1,l2,l3])

            ID2 = len(get_curvilinear_coordinate(test_cl))
            idmid = int( (ID1 + ID2)/2.)

            for p in range(voronoi_clipped.GetNumberOfPoints()):
                cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

                if cl_id < idmid:
                    dist = dx_p1 * (idmid**2-cl_id**2) / (idmid**2-ID1**2)
                else:
                    if cl_id <= (ID2 - 4):
                        dist = -dx_p1 * (cl_id-idmid)**(0.5) / (ID2-idmid)**(0.5)
                    else: 
                        # Move opthalmic artery based on its discovery ID
                        if clip_ID < idmid_0:
                            cl_id =  clip_ID - ID1_0  
                            dist = dx_p1 * (idmid**2-cl_id**2) / (idmid**2-ID1**2)
                        else:
                            cl_id = clip_ID - ID1_0 
                            dist = -dx_p1 * (cl_id-idmid)**(0.5) / (ID2-idmid)**(0.5)

                points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
                verts.InsertNextCell(1)
                verts.InsertCellPoint(p)

        else:
            ID2 = len(get_curvilinear_coordinate(patch_cl))
            idmid = int( (ID1 + ID2)/2.)
            ID1=ID2=idmid
            for p in range(voronoi_clipped.GetNumberOfPoints()):
                cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))
                if cl_id < ID1:
                    ID1 = cl_id

                if cl_id > ID2:
                    ID2 = cl_id
    

            for p in range(voronoi_clipped.GetNumberOfPoints()):
                cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

                if cl_id < idmid:
                    dist = dx_p1 * (idmid**2-cl_id**2) / (idmid**2-ID1**2)
                else:
                    dist = -dx_p1 * (cl_id-idmid)**(0.5) / (ID2-idmid)**(0.5)

                points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
                verts.InsertNextCell(1)
                verts.InsertCellPoint(p)

    else:
        for p in range(voronoi_clipped.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

            if cl_id <= ID1:
                dist = dx_p1
            else:
                dist =  -dx_p1
            
            points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newDataSet.SetPoints(points)
    newDataSet.SetVerts(verts)
    newDataSet.GetPointData().AddArray(voronoi_clipped.GetPointData().GetArray(radiusArrayName))

    return newDataSet


def move_voro(clip_voro, old_cl, ID1_0, clip_ID, dx,eye=False):

    centerline_loc = get_locator(old_cl)

    newDataSet = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()
    if eye == True:
        l1 = ExtractSingleLine(old_cl, 0)
        l2 = ExtractSingleLine(old_cl, 1)
        l3 = ExtractSingleLine(old_cl, 2)
        test_cl = merge_data([l1,l2,l3])

        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(test_cl))
        IDmid = int( (ID1 + ID2)/2.)
        I1 = I2 = IDmid

        for p in range(clip_voro.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(clip_voro.GetPoint(p))
            if cl_id < I1:
                I1 = cl_id

            if cl_id > I2 and cl_id <= (ID2-4):
                I2 = cl_id

        for p in range(clip_voro.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(clip_voro.GetPoint(p))

            if cl_id <= (ID2-4):
                dist = 4 * dx * (cl_id - I1)*(I2 - cl_id) / ( I2 - I1)**2
            else:
                cl_id = clip_ID - ID1_0 
                dist = 4 * dx * (cl_id - ID1)*(ID2 - cl_id) / ( ID2 - ID1)**2

            points.InsertNextPoint(np.asarray(clip_voro.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    else:
        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(old_cl)) - 2
        IDmid = int( (ID1 + ID2)/2.)
        ID1 = ID2 = IDmid

        for p in range(clip_voro.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(clip_voro.GetPoint(p))
            if cl_id < ID1:
                ID1 = cl_id

            if cl_id > ID2:
                ID2 = cl_id
    
        for p in range(clip_voro.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(clip_voro.GetPoint(p))

            dist = 4 * dx * (cl_id - ID1)*(ID2 - cl_id) / ( ID2 - ID1)**2

            points.InsertNextPoint(np.asarray(clip_voro.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newDataSet.SetPoints(points)
    newDataSet.SetVerts(verts)
    newDataSet.GetPointData().AddArray(clip_voro.GetPointData().GetArray(radiusArrayName))
    newVoronoi = newDataSet
    return newVoronoi


if  __name__ == "__main__":
    smooth, basedir, case, alpha, beta = read_command_line()
    folders = sorted([folder for folder in listdir(basedir) if folder[:2] in ["P0"]])
    name = "surface"
    growDir = "both"
    rad = 0.2
    filepath = "plane_curv/alphabeta_values.txt"
    #filepath = "alphabeta_values_neigb_%s_radius_%s.txt" % (NUM,rad)
    myfile = open(filepath,"r")


    for folder in folders[3:]:
        info = myfile.readline().split()
        info = myfile.readline().split()
        info = myfile.readline().split()
        info = myfile.readline().split()
        info = myfile.readline().split()
        info = myfile.readline().split()
        for i in range(2):
            info = myfile.readline().split()
            alpha = float(info[-2].split("=")[1])
            beta = float(info[-1].split("=")[1])
            #alpha = 0.1
            #beta = 0.5
            print alpha,beta
            
            print("==== Working on case %s ====" % folder)
            case = path.join(basedir,folder)
            #check = path.join(case, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
            #if not path.isfile(check):
            ##if folder =="P0157":
            reconst_and_move_voro(case ,smooth, name, alpha, beta)
        break




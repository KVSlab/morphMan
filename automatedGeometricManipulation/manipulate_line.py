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
from moveandmanipulatetools import *

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


def reconst_and_move_voro(dirpath, smooth, name, alpha=0.0, beta=0.0):
    case = dirpath[-5:]
    # Input filenames
    model_path = path.join(dirpath, name, "model.vtp")

    # Output names
    model_smoothed_path = path.join(dirpath, name, "model_smoothed.vtp")
    model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
    model_new_surface_tmp = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s_tmp.vtp" % (alpha, beta))
    new_centerlines_path = path.join("./modelsangle", "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    model_new_surface = path.join("./modelsangle", "%s_alpha_%s_beta_%s.vtp" % (case, alpha,beta))

    new_centerlines_path = path.join("./testmodels", "%s_centerlines_alpha_%s_beta_%s_.vtp" % (case, alpha, beta))
    model_new_surface = path.join("./testmodels", "%s_alpha_%s_beta_%s_.vtp" % (case, alpha,beta))
    new_centerlines_path = path.join("./animodels", "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    new_centerlines_path_tmp = path.join("./animodels", "%s_centerlines_alpha_%s_beta_%s_tmp.vtp" % (case, alpha, beta))
    model_new_surface = path.join("./animodels", "%s_alpha_%s_beta_%s.vtp" % (case, alpha,beta))

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

    #voronoi = voronoi if not smooth else ReadPolyData(voronoi_smoothed_path)
    voronoi = SmoothClippedVoronoiDiagram(voronoi, centerlines_complete, 0.25)
    #WritePolyData(voronoi_smoothed, voronoi_smoothed_path)
    #WritePolyData(voronoi, "vorotmp.vtp")

    # Set clipping points
    div_points = np.asarray(clipping_points)
    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)
    clip_points = points


    centerlines_in_order = sort_centerlines(centerlines_complete)


    # Check if case includs the opthalmic artery
    eye = False
    eye, clip_ID, centerlines_in_order, eyeline = find_ophthalmic_artery(centerlines_in_order, clipping_points)

    if eye == True:
        print("Clipping opthamlic artery")
        patch_eye = clip_eyeline(eyeline, clipping_points[0], clip_ID)

    # Get location of clipping points
    line = ExtractSingleLine(centerlines_in_order, 0)
    locator = get_locator(line)
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)

    # Clip centerline
    print("Clipping centerlines.")
    patch_cl = CreateParentArteryPatches(centerlines_in_order, clip_points)
    clipped_curve = ExtractSingleLine(centerlines_in_order, 0, startID=ID1, endID=ID2)
    WritePolyData(patch_cl, centerline_clipped_path)
    
    if eye == True:
        eyeline_end = ExtractSingleLine(patch_eye,1)
        clipped_curve = merge_data([clipped_curve, eyeline_end])

    WritePolyData(clipped_curve, centerline_clipped_part_path )


    # Clipp the voronoi diagram
    # And keep the clipped part
    print("Computing clipped voronoi diagrams")
    #if not path.exists(voronoi_clipped_path):
    masked_voronoi_clip  = MaskVoronoiDiagram(voronoi, clipped_curve)
    voronoi_clipped_part = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi_clip)
    WritePolyData(voronoi_clipped_part, voronoi_clipped_part_path)
    #else:
    #    voronoi_clipped_part = ReadPolyData(voronoi_clipped_part_path)

    #if not path.exists(voronoi_clipped_part_path):
    masked_voronoi  = MaskVoronoiDiagram(voronoi, patch_cl)
    voronoi_clipped = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
    WritePolyData(voronoi_clipped, voronoi_clipped_path)
    #else:
    #    voronoi_clipped = ReadPolyData(voronoi_clipped_path)
    
    # Extract translation vectors
    print("Computing translation directions.")
    growDir = "horizont"
    middle_points, middleIds = get_spline_points(line, beta, growDir, case, clip_points,voronoi)

    # Move anterior bend horizontally 

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
        WritePolyData(voronoi_clipped, "a_voronoi_clipped.vtp")
        WritePolyData(voronoi_clipped_part, "a_voronoi_clipped_part.vtp")


        # Move line manually for postprocessing
        # Update inlet and outlet positions
        print("Adjusting line manually")
        line1 = ExtractSingleLine(patch_cl, 0)
        lines = []

        n = patch_cl.GetNumberOfCells()
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
        # Move line manually for postprocessing
        # Update inlet and outlet positions
        print("Adjusting line manually")
        line1 = ExtractSingleLine(patch_cl, 0)
        lines = []

        n = patch_cl.GetNumberOfCells()
        for i in range(n-1):
            lines.append(ExtractSingleLine(patch_cl, i+1))
            
        patch_cl_start = line1
        patch_cl_ends = []
        for line_n in lines:
            patch_cl_end = line_n
            newpoints.append(patch_cl_end.GetPoint(patch_cl_end.GetNumberOfCells() - 1))
            patch_cl_ends.append(patch_cl_end)


        clipped_part_new = clipped_curve

        if eye:
            newpoints.append(clipped_part_new.GetPoint(clipped_part_new.GetNumberOfCells() - 1))
        
        newpoints.append(patch_cl_start.GetPoint(0))

        newVoronoi = voronoi
        new_surface = surface

    print("Writing new surface (tmp)")
    WritePolyData(new_surface, model_new_surface_tmp)

    if alpha != 0.0:
        if beta == 0.0:
            print("No horizontal movement. Initiating vertical movement")
            new_centerline = centerlines_complete
            WritePolyData(new_centerline, new_centerlines_path) 
        else:
            print "Creating new_centerline_complete.vtp of horizontally moved model"
            new_centerline = makeCenterline(model_new_surface_tmp, new_centerlines_path, smooth=False, resampling=False, newpoints=newpoints, recompute=True)
            WritePolyData(new_centerline, new_centerlines_path_tmp) 

        # Vertical movement
        print("Moving geometry vertically")
        translate_vertically(new_surface, new_centerline, clip_points, alpha , dirpath, newVoronoi, voronoi_clipped, voronoi_clipped_part, eye, newpoints, case)
    
    else:
        # Write a new surface from the new voronoi diagram
        WritePolyData(new_surface, model_new_surface)
    return new_surface



def translate_vertically(surface, centerline, clip_points, alpha, dirpath, voronoi, voronoi_clipped, voronoi_clipped_part, eye, oldpoints, case):
    # Filenames
    #new_centerlines_path = path.join(dirpath, name, "new_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta))
    #model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
    new_centerlines_path = path.join("./modelsangle", "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    model_new_surface = path.join("./modelsangle", "%s_alpha_%s_beta_%s.vtp" % (case, alpha,beta))

    new_centerlines_path = path.join("./animodels", "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    model_new_surface = path.join("./animodels", "%s_alpha_%s_beta_%s.vtp" % (case, alpha,beta))

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
    locator = get_locator(ExtractSingleLine(centerlines_in_order,0))
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)
    clipped_curve = ExtractSingleLine(centerlines_in_order, 0, startID=ID1, endID=ID2)

    #clipped_curve = get_clipped_curve(centerlines_in_order, clipping_points) 

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
    WritePolyData(moved_clip, "voromovedpartclip_1.vtp")
    WritePolyData(voronoi_clipped, "voromovedclip_1.vtp")

    # Move centerline manually for postprocessing  
    print("Adjusting line manually")
    if eye:
        newpoints = oldpoints[:-2] + [oldpoints[-1]]
    else:
        newpoints = oldpoints

    clipped_part_new = move_dots_vertical(clipped_curve, dx, ID1, clip_ID, eye)

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
            # Find limits of anterior bend curve
            l1 = ExtractSingleLine(patch_cl, 0)
            ID2 = len(get_curvilinear_coordinate(l1))
            idmid = int( (ID1 + ID2)/2.)

            for p in range(voronoi_clipped.GetNumberOfPoints()):
                cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

                if cl_id < idmid:
                    dist = dx_p1 * (idmid**2-cl_id**2) / (idmid**2-ID1**2)
                else:
                    if cl_id <= (ID2-1):
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
        for p in range(voronoi_clipped.GetNumberOfPoints()+1):
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

        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(l1)) 
        IDmid = int( (ID1 + ID2)/2.)
        I1 = I2 = IDmid

        for p in range(clip_voro.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(clip_voro.GetPoint(p))
            if cl_id < I1:
                I1 = cl_id

            if cl_id > I2 and cl_id <= (ID2-1):
                I2 = cl_id

        for p in range(clip_voro.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(clip_voro.GetPoint(p))

            if cl_id <= (ID2-1):
                dist = 4 * dx * (cl_id - I1)*(I2 - cl_id) / ( I2 - I1)**2
            else:
                cl_id = clip_ID - ID1_0 #+10#+ int((ID2 - (clip_ID - ID1_0) )*0.4)
                dist = 4 * dx * (cl_id - ID1)*(ID2 - cl_id) / ( ID2 - ID1)**2
                #cl_id = clip_ID - ID1_0 
                #dist = 4 * dx * (cl_id - ID1)*(ID2 - cl_id) / ( ID2 - ID1)**2

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
    #filepath = "plots_curv/optimal_alphabeta.txt"
    filepath = "plots_angle/optimal_alphabeta.txt"
    #filepath = "alphabeta_values_neigb_%s_radius_%s.txt" % (NUM,rad)


    b1 =  0.5
    a1 = 0.4
    # P0220
    #a1=0.145484949833
    #b1=0.200836120401
    #a2=0.0290969899666
    #b2=-0.164882943144

    # P0250
    a1 = 0.0451505016722
    b1 =0.145150501672
    a2 =-0.123411371237
    b2 =-0.0875585284281
    #a2 =0.101337792642
    #b2 =0.113545150502
    #a1=-0.0872909698997 
    #b1=-0.145652173913
    # P0207
    a1= 0.117391304348
    b1=-0.0961538461538
    a2=0.129431438127
    b2=0.0769230769231

    # P=134
    #a1=0.0973244147157
    #b1=0.115384615385
    #a2=0.145484949833
    #b2=-0.0384615384615

    # NOTE: ANGLE P0086
    a1=-0.2 
    b1=0.4506093489148581
    a2=0.844908180300501 
    b2=-0.0009933222036727918

    a1 = 0.543739
    b1 = -0.175112
    a2=0.3899142
    b2 = 0.41139192
    
    a_ = [a1,a2]
    b_ = [b1,b2]
    a2 = [0.01, 0.1, 0.2, 0.3, 0.4, 0.5 ,0.6 ,0.7 ,0.8, 0,0,0,0,0,0,0,0, 0.1] 

    #a1=[0.05, 0.1,0.15 ,
    a1=[0.2,0.25,0.3,0.35,0.4, 0.45,0.5,0.55]
    #b1=[0.1, 0.2 ,0.3 ,
    b1=[0.4 ,0.5 ,0.6,0.7 , 0.8,0.9,1.0, 1.1]
    b2 = [0   ,0   ,0   ,0  ,0   ,0   ,0    ,0    ,0  , 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.1]
    a_ = a1+a2
    b_ = b1+b2

    a1=[0.05, 0.1,0.15, 0.2,0.25,0.3,0.35,0.4, 0.45,0.5,0.55]
    b1=[0.1, 0.2 ,0.3 , 0.4 ,0.5 ,0.6,0.7 , 0.8,0.9,1.0, 1.1]
    b2 = [0   ,0   ,0   ,0  ,0   ,0   ,0    ,0    ,0  , 0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.1]
    a_ = a1+a2
    b_ = b1+b2
    print a_
    print b_



    #b_= [0.1, 0.2,0.3,0.4,0.5,      0.6,0.7,0.8,0.9, 1, 1.1, 1.2, 1.3]
    #a_= [-0.1,-0.2, -0.1, 0.0, 0.1, 0.2,0.3,0.4,0.5,0.4,0.3, 0.2, 0.1]
    print len(a_), len(b_)
    myfile = open(filepath,"r")
    info = myfile.readline().split()
    info = myfile.readline().split()
    info = myfile.readline().split()
    info = myfile.readline().split()
    info = myfile.readline().split()
    info = myfile.readline().split()
    for folder in folders[0:1]: 
    #for folder in [folders[2], folders[6]]:
        #for alpha, beta in zip(a_,b_):

        #info = myfile.readline().split()
        #info = myfile.readline().split()
        #info = myfile.readline().split()
        #info = myfile.readline().split()
        #info = myfile.readline().split()
        k=0
        for i in range(len(a_)):
            info = myfile.readline().split()
            #alpha = float(info[-2].split("=")[1])
            #beta = float(info[-1].split("=")[1])
            #alpha = 0.543739
            #beta = -0.175
            alpha = a_[i]
            beta = b_[i]
            print alpha,beta
            
            print("==== Working on case %s ====" % folder)
            case = path.join(basedir,folder)
            #check = path.join(case, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
            #if not path.isfile(check):
            ##if folder =="P0157":
            reconst_and_move_voro(case ,smooth, name, alpha, beta)




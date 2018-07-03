from common import *
from argparse import ArgumentParser
from os import path, remove, listdir
from subprocess import STDOUT, check_output
from IPython import embed
from math import isnan

import matlab.engine
import operator
import sys
import math
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema, resample
from scipy.ndimage.filters import gaussian_filter as gauss
from mpl_toolkits.mplot3d import Axes3D

# Local import
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from matplotlib import rc, rcParams

rc('text', usetex=True)
rc('font', family='serif')
rcParams['axes.linewidth'] = 2.0
rcParams['figure.figsize'] = 15, 10

# Arrays for storing angle variation
a1=[];A1=[]
a2=[];A2=[]
a3=[];A3=[]
a4=[];A4=[]
a5=[];A5=[]
m1=[];M1=[]
m2=[];M2=[]
m3=[];M3=[]
m4=[];M4=[]
s1=[];s2=[];s3=[];s4=[];
S1=[];S2=[];S3=[];S4=[]
dp1=[];dp2=[];dp3=[];dp4=[];dp5=[]
dm1=[];dm2=[];dm3=[];dm4=[]
ds1=[];ds2=[];ds3=[]; ds4=[]
svd_d0 = []; svd_d1 = []
z1 = []; z2 = [];z3=[];z4=[]
Z1 = []; Z2 = [];Z3=[];Z4=[]
dz1 = []; dz2=  [];dz3=[];dz4=[]
old_deg = []
new_deg = []

old_deg1 = []
old_deg2 = []
old_deg3 = []
new_deg1 = []
new_deg2 = []
new_deg3 = []


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('-d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('--case', type=str, default=None, help="Choose case")
    parser.add_argument('-s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram (surface) should be" + \
                        "smoothed before it is manipulated", metavar="smooth")

    parser.add_argument('-k', '--curvature', type=bool, default=False,
                        help="Compute curvature variation", metavar="curvature")
    parser.add_argument('-t', '--angle', type=bool, default=False,
                        help="Compute angle variation", metavar="angle")
    parser.add_argument("-a", "--alpha", type=float, default=None,
                        help="Compression factor in vertical direction, ranging from -1.0 to 1.0")
    parser.add_argument("-b", "--beta", type=float, default=None,
                        help="Compression factor in horizontal direction, ranging from -1.0 to 1.0")

    args = parser.parse_args()

    return args.smooth, args.dir_path, args.case, args.curvature, args.angle, args.alpha, args.beta


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
    particles = path.join(dirpath, "optimal_bog.particles")
    all_points = np.loadtxt(particles)
    clipping_points = all_points[2:] 
    return clipping_points


def compute_angle(dirpath, smooth, name, alpha,beta, proj=False):
    case = dirpath[-5:]
    # Input filenames
    model_path = path.join(dirpath, name, "model.vtp")

    # Output names
    model_smoothed_path = path.join(dirpath, name, "model_smoothed.vtp")
    model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha, beta))
    new_centerlines_path = path.join(dirpath, name, "new_centerlines_complete_alpha_%s_beta_%s.vtp" % (alpha, beta))
    dirpath2 = "modelstable"
    model_new_surface = path.join(dirpath2, "%s_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    new_centerlines_path = path.join(dirpath2, "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))

    # Centerliens
    centerline_complete_path = path.join(dirpath, name, "centerline_complete.vtp")
    centerline_clipped_path = path.join(dirpath, name,  "centerline_clipped.vtp")
    centerline_clipped_part_path = path.join(dirpath, name,  "centerline_clipped_part.vtp")
    centerline_clipped_part1_path = path.join(dirpath, name,  "centerline_clipped_end_part.vtp")
    centerline_new_path = path.join(dirpath, name, "centerline_interpolated.vtp")
    centerlines_ordered_path = path.join(dirpath, name, "centerline_ordered.vtp")

    # Find endID from landmarking
    siphon_path = path.join(dirpath, "surface", "carotid_siphon.vtp")

    # Voronoi diagrams
    voronoi_path = path.join(dirpath, name, "model_voronoi.vtp")
    voronoi_smoothed_path = path.join(dirpath, name,  "model_voronoi_smoothed.vtp")
    voronoi_clipped_path = path.join(dirpath, name, "model_voronoi_clipped.vtp")
    voronoi_clipped_part_path = path.join(dirpath, name, "model_voronoi_clipped_part.vtp")
   
    # Extract Clipping points
    compute_old = False
    clipping_points = get_clipping_points(dirpath)
    div_points = np.asarray(clipping_points)
    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)
    clip_points = points

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
    
    # Compute centerline
    centerlines_complete = compute_centerlines(inlet, outlets,
                                               centerline_complete_path,
                                               capped_surface, resampling=0.1)

    voronoi = makeVoronoiDiagram(surface, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerlines_complete, 0.25)
        WritePolyData(voronoi_smoothed, voronoi_smoothed_path)
        surface_smoothed = create_new_surface(voronoi_smoothed)
        WritePolyData(surface_smoothed, model_smoothed_path)

    voronoi = voronoi if not smooth else ReadPolyData(voronoi_smoothed_path)


    # Special cases including the opthalmic artery
    eye = False
    eye, clip_ID, centerlines_complete, eyeline = find_ophthalmic_artery(centerlines_complete, clipping_points)

    # FInd longest line:
    lines = []
    n = centerlines_complete.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines_complete, i))

    longest = [lines[0]]
    lenlong = get_curvilinear_coordinate(longest[0])
    for i in range(1,n):
        tmplong = get_curvilinear_coordinate(lines[i])
        if len(tmplong) > len(lenlong):
            lenlong = tmplong
            longest.insert(0, lines[i])
        else: 
            longest.append(lines[i])


    # Move new line manually
    centerlines_in_order = merge_data(longest)
    line = ExtractSingleLine(centerlines_in_order, 0)
    patch_cl = CreateParentArteryPatches(centerlines_in_order, clip_points)
    print("Getting clipped curve.")
    clipped_curve = get_clipped_curve(centerlines_in_order, clipping_points) 
    
    patch_start = ExtractSingleLine(patch_cl, 0)
    patch_ends = []
    for i in range(1,n+1):
        patch_ends.append(ExtractSingleLine(patch_cl, i))

    # Find ID of middle pooint:
    growDir_tmp = "horizont"
    middle_points, middleIds = get_spline_points(line, beta, growDir_tmp, dirpath, clip_points,voronoi)

    locator = get_locator(line)
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)
    ID_mid = int((ID1 + ID2) / 2.)
    dx_p1 = middle_points[0] - p1
    dx_p2 = middle_points[-1] - p2

    # Compute new centerline by manual displacement
    print("Moving centerline manually")
    patchline_1 = patch_start
    patchline_2 = patch_ends[0]
    patch_cl_new1 = move_line_horo(patchline_1, ID1, ID2, dx_p1, clip=False, side="left")
    patch_cl_new2 = move_line_horo(patchline_2, ID1, ID2, dx_p1, clip=False, side="right")
    clipped_part_new = move_line_horo(clipped_curve, ID1, ID2, dx_p1, clip=True, eye=eye)
    new_centerline = merge_data([patch_cl_new1, clipped_part_new, patch_cl_new2])
    clipped_part_new = connect_line(clipped_part_new)
    growDir = "vertical"
    middle_points, middleIds, dx = get_spline_points(new_centerline, alpha, growDir, dirpath, clip_points,voronoi)
    clipped_part_new = move_dots_vertical(clipped_part_new, dx, eye=eye)
    new_centerline = merge_data([patch_cl_new1, clipped_part_new, patch_cl_new2])
    new_centerline = connect_line(new_centerline)


    print("Spline centerline and add curvature array") 
    siphon = ReadPolyData(siphon_path)
    # Spline line as in Landmarking
    nknots = 11
    line_sp, curv_sp = spline_line(siphon, get_curv=True, isline=True, nknots=nknots)
    curv_sp= resample(curv_sp, line_sp.GetNumberOfPoints())
    
    # Smooth line with discrete derivatives
    neigh = 30
    
    # Spline moved line 
    locator = get_locator(new_centerline)
    endID = siphon.GetNumberOfPoints() - 1
    endID_new  = locator.FindClosestPoint(line_sp.GetPoint(endID))
    newsiphon = cutLine(new_centerline, endID_new)
    newline_sp, newcurv_sp = spline_line(newsiphon, get_curv=True, isline=True, nknots=nknots)
    
    # Smooth line with discrete derivatives
    #line_d, curv_d =  discrete_geometry(dirpath, siphon, neigh=neigh)
    #newline_d, newcurv_d =  discrete_geometry(dirpath, newsiphon, neigh=neigh)

    # Extract anterior bend 
    #line_sp = line
    locator = get_locator(line_sp)
    ID1  = locator.FindClosestPoint(p1)
    ID2  = locator.FindClosestPoint(p2)

    locator = get_locator(newline_sp)
    newID1  = locator.FindClosestPoint(p1)
    newID2  = locator.FindClosestPoint(p2)
    old_siphon = line_sp
    new_siphon = newline_sp
    siphon = ExtractSingleLine(old_siphon, 0, startID=ID1, endID=ID2)
    moved_siphon = ExtractSingleLine(new_siphon, 0, startID=newID1, endID=newID2)
    #WritePolyData(siphon, "angledata/siphon_before_%s.vtp" % case)
    #WritePolyData(moved_siphon, "angledata/siphon_after_%s.vtp" % case)

    new_p1 = newline_sp.GetPoints().GetPoint(newID1)
    new_p2 = newline_sp.GetPoints().GetPoint(newID2)
    """

    # Map MISR values to old and new splined anterior bend
    anterior_bend =ExtractSingleLine( line, 0, startID=ID1, endID=ID2)
    M = anterior_bend.GetNumberOfPoints()
    M1 = moved_siphon.GetNumberOfPoints()
    misrArray = get_vtk_array(radiusArrayName, 1 ,M)
    newmisrArray = get_vtk_array(radiusArrayName, 1 ,M1)
    MISR = []
    for i in range(M):
        misr = anterior_bend.GetPointData().GetArray(radiusArrayName).GetTuple(i)
        MISR.append(misr[0])
        misrArray.SetTuple(i, misr)

    MISR = resample(MISR, M1)
    for i in range(M1):
        newmisrArray.SetTuple(i, (MISR[i],))

    siphon.GetPointData().AddArray(misrArray)
    moved_siphon.GetPointData().AddArray(newmisrArray)
    """

    """ Find midpoint / Max curv point / Fractions for angle measuremetn"""
    methods = ["odrline", "MISR", "frac", "plane", "itplane", "itplane_clip", "maxcurv", "smooth", "discrete", "maxdist"] 
    methods = ["plane"]#
    
    cutcurv = curv_sp[ID1:ID2]
    newcutcurv = newcurv_sp[newID1:newID2]

    #cutcurv_d = curv_d[ID1:ID2]
    #newcutcurv_d = newcurv_d[newID1:newID2]

    import operator
    if proj:
        print "Computing 2D Angles"
    else:
        print "Computing 3D Angles"
        

    for i,method in enumerate(methods):
        if method == "plane":
            maxP, maxID     = find_furthest_points(dx, p1,p2, siphon, show=False)
            newmaxP, newmaxID = find_furthest_points(dx,new_p1,new_p2, moved_siphon, show=False)

        elif method in ["itplane", "itplane_clip"]:
            maxP, maxID     = find_furthest_points(dx, p1,p2, siphon, show=False)
            newmaxP, newmaxID = find_furthest_points(dx,new_p1,new_p2, moved_siphon, show=False)

            T1 = get_array("FrenetTangent", siphon, k=3)
            T2 = get_array("FrenetTangent", moved_siphon, k=3)
            
            p1_1, p1_id = find_closest_point(T1[-1],0, maxID, p2, siphon, show=False)
            p2_2, p2_id = find_closest_point(T1[0],maxID, siphon.GetNumberOfPoints(), p1, siphon, show=False)

            newp1_1, np1_id = find_closest_point(T2[-1],0, newmaxID , new_p2, moved_siphon, show=False)
            newp2_2, np2_id = find_closest_point(T2[0],newmaxID, moved_siphon.GetNumberOfPoints(), new_p1, moved_siphon, show=False)

            N1 = get_array("FrenetBinormal", siphon, k=3)[p1_id]
            N2 = get_array("FrenetBinormal", moved_siphon, k=3)[np1_id]

            dP = p1_1 - p2_2
            dnewP = newp1_1 - newp2_2 

            normal = np.cross(dP, N1)
            newnormal = np.cross(dnewP, N2)
            
            maxP, maxID     = find_furthest_points(normal, p1_1,p2_2, siphon, show=False)
            newmaxP, newmaxID = find_furthest_points(newnormal,newp1_1,newp2_2, moved_siphon, show=False)

        elif method == "maxcurv":
            maxID,v = max(enumerate(cutcurv), key=operator.itemgetter(1))
            newmaxID,v = max(enumerate(newcutcurv), key=operator.itemgetter(1))
            maxP = line_sp.GetPoint(ID1 + maxID)
            newmaxP = newline_sp.GetPoint(newID1+newmaxID)

        elif method == "smooth":
            allmaxcurv = argrelextrema(cutcurv, np.greater)[0]
            allnewmaxcurv = argrelextrema(newcurv_sp, np.greater)[0]

            tmpcurv = cutcurv
            while len(allmaxcurv) > 2:
                tmpcurv = gauss(tmpcurv, 2)
                allmaxcurv = argrelextrema(tmpcurv, np.greater)[0]

            tmpnewcurv = newcutcurv
            while len(allnewmaxcurv) > 2:
                tmpnewcurv = gauss(tmpnewcurv, 2)
                allnewmaxcurv = argrelextrema(tmpnewcurv, np.greater)[0]

            maxID = allmaxcurv[0]
            newmaxID = allnewmaxcurv[0]

        elif method == "discrete":
            maxID,v = max(enumerate(cutcurv_d), key=operator.itemgetter(1))
            newmaxID,v = max(enumerate(newcutcurv_d), key=operator.itemgetter(1))


        elif method == "maxdist":
            normP1 = [la.norm(np.array(p1) - np.array(siphon.GetPoint(i))) for i in range(siphon.GetNumberOfPoints())]
            normP2 = [la.norm(np.array(p2) - np.array(siphon.GetPoint(i))) for i in range(siphon.GetNumberOfPoints()-1,-1,-1)]
            maxID=0;  max_dist = 0
            for i,n1 in enumerate(normP1):
                for j,n2 in enumerate(normP2):
                    dist = n1**2 + n2**2
                    if dist > max_dist:
                        max_dist = dist; maxID=i
            
            newnormP1 = [la.norm(np.array(new_p1) - np.array(moved_siphon.GetPoint(i))) for i in range(moved_siphon.GetNumberOfPoints())]
            newnormP2 = [la.norm(np.array(new_p2) - np.array(moved_siphon.GetPoint(i))) for i in range(moved_siphon.GetNumberOfPoints()-1,-1,-1)]
            newmaxID= 0; new_max_dist = 0
            for i,n1 in enumerate(newnormP1):
                for j,n2 in enumerate(newnormP2):
                    dist = n1**2 + n2**2
                    if dist > new_max_dist:
                        new_max_dist = dist; newmaxID=i


        #print("Computing angles")
        if method == "odrline":
            limits = ["cumulative" , "sd"] 
            pA = pB = newpA = newpB = np.zeros(3)
            for limit in limits:
                p01, d1, p02, d2, curvlineold = odr_line(ID1, ID2, line_sp, curv_sp, limit)
                p01, newd1, p02, newd2, curvlinenew = odr_line(newID1, newID2, newline_sp, newcurv_sp, limit)
                WritePolyData(curvlineold, "angledata/%s_curlineold.vtp" % case)
                WritePolyData(curvlinenew, "angledata/%s_curlinenew.vtp" % case)

                deg,l1,l2 = find_angle_odr(d1,d2,proj)
                newdeg,nl1,nl2 = find_angle_odr(newd1,newd2,proj)

                save_degree(deg, pA,pB, p1,p2,l1,l2, method, limit, case, "before",proj)
                save_degree(newdeg, newpA,newpB, new_p1, new_p2, nl1,nl2, method, limit,case, "after", proj)

                old_deg.append(deg)
                new_deg.append(newdeg)

        elif method == "MISR":
            multiples = [1,1.5,2,2.5]
            for param in multiples:
                N1 = siphon.GetNumberOfPoints()
                N2 = moved_siphon.GetNumberOfPoints()
                rad1 = siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(0)
                rad2 = siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(N1-1)
                newrad1 = moved_siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(0)
                newrad2 = moved_siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(N2-1)

                pA, rA = move_past_sphere(siphon, p1, rad1, 0, step=1, stop=N1-1, X=param)
                pB, rB = move_past_sphere(siphon, p2, rad2, N1-1, step=-1, stop=0, X=param)
                newpA, rA = move_past_sphere(moved_siphon, new_p1, newrad1, 0, step=1, stop=N2-1, X=param)
                newpB, rB = move_past_sphere(moved_siphon, new_p2, newrad2, N2-1, step=-1, stop=0, X=param)

                deg,l1,l2 = find_angle(pA, pB, p1,p2, proj)
                newdeg,nl1,nl2 = find_angle(newpA, newpB, new_p1,new_p2, proj)
                
                save_degree(deg, pA,pB, p1,p2, l1,l2, method, param, case, "before",proj)
                save_degree(newdeg, newpA,newpB, new_p1, new_p2, nl1,nl2, method, param,case, "after",proj)

                old_deg.append(deg)
                new_deg.append(newdeg)

        else: 
            if method == "frac":
                n_values = [5]#[3,5,5,7,7]
                l = [2]#[1,2,1,1,2]
                r = [3]##[2,3,4,6,5]
                for i in range(len(n_values)):
                    dX = 1. / n_values[i]
                    frac = "%sdiv%s"% (l[i] ,n_values[i])
                    IDA = int(ID1 + (ID2 - ID1) * l[i] * dX) 
                    IDB = int(ID1 + (ID2 - ID1) * r[i] * dX) 
                    pA = line_sp.GetPoints().GetPoint(IDA)
                    pB = line_sp.GetPoints().GetPoint(IDB)

                    IDA = int(newID1 + (newID2 - newID1) * l[i] * dX) 
                    IDB = int(newID1 + (newID2 - newID1) * r[i] * dX) 
                    newpA = newline_sp.GetPoints().GetPoint(IDA)
                    newpB = newline_sp.GetPoints().GetPoint(IDB)

                    #deg,l1,l2 = find_angle(pA, pB, p1,p2, proj)
                    newdeg,nl1,nl2 = find_angle(newpA, newpB, new_p1,new_p2, proj)
                    #FIXME 
                    return newdeg
                    
                    save_degree(deg, pA,pB, p1,p2, l1,l2, method, frac, case, "before",proj)
                    save_degree(newdeg, newpA,newpB, new_p1, new_p2, nl1,nl2, method, frac,case, "after",proj)

                    old_deg.append(deg)
                    new_deg.append(newdeg)

            elif method in ["plane", "itplane", "itplane_clip", "maxcurv", "smooth", "discrete", "maxdist"]: 
                frac_vals = [ 2. / 5. ]
                for frac in frac_vals:

                    if method == "itplane_clip":
                        IDmid = (p2_id - p1_id)/2.#siphon.GetNumberOfPoints()
                        newIDmid = (np2_id - np1_id)/2.#moved_siphon.GetNumberOfPoints()
                        if maxID > IDmid: 
                            IDA = int((maxID-p1_id)*frac)
                            IDB = int((maxID-p1_id)*(1 + (1-frac)))
                            pA = siphon.GetPoints().GetPoint(IDA+p1_id)
                            pB = siphon.GetPoints().GetPoint(IDB+p1_id)
                        else:
                            IDB = int((p2_id - maxID)*(1 + (1-frac)))
                            IDA = int((p2_id - maxID)*frac)
                            pA = siphon.GetPoints().GetPoint(IDA)
                            pB = siphon.GetPoints().GetPoint(IDB)

                        if newmaxID > newIDmid:
                            IDA = int((newmaxID-np1_id)*frac)
                            IDB = int((newmaxID-np1_id)*(1 + (1-frac)))
                            newpA = moved_siphon.GetPoints().GetPoint(IDA+np1_id)
                            newpB = moved_siphon.GetPoints().GetPoint(IDB + np1_id)
                        else:
                            IDA = int((np2_id - newmaxID)*frac)
                            IDB = int((np2_id - newmaxID)*(1 + (1-frac)))
                            newpA = moved_siphon.GetPoints().GetPoint(IDA)
                            newpB = moved_siphon.GetPoints().GetPoint(IDB)


                        deg,l1,l2 = find_angle(pA, pB, p1_1,p2_2, proj)
                        newdeg,nl1,nl2 = find_angle(newpA, newpB, newp1_1,newp2_2, proj)
                        save_degree(deg, pA,pB, p1_1,p2_2, l1,l2, method, frac, case, "before",proj)
                        save_degree(newdeg, newpA,newpB, newp1_1, newp2_2, nl1,nl2, method, frac,case, "after",proj)

                        old_deg.append(deg)
                        new_deg.append(newdeg)
                    else:
                        IDA = int(maxID*frac)
                        IDB = int(maxID*(1 + (1-frac)))
                        pA = siphon.GetPoints().GetPoint(IDA)
                        #PM = siphon.GetPoints().GetPoint(int(maxID))
                        pB = siphon.GetPoints().GetPoint(IDB)

                        deg,l1,l2 = find_angle(pA, pB, p1,p2, proj)
                        save_degree(deg, pA,pB, p1,p2, l1,l2, method, frac, case, "before",proj)
                        embed()

                        IDA = int(newmaxID*frac)
                        IDB = int(newmaxID*(1 + (1-frac)))
                        newpA = moved_siphon.GetPoints().GetPoint(IDA)
                        newpB = moved_siphon.GetPoints().GetPoint(IDB)

                        newdeg,nl1,nl2 = find_angle(newpA, newpB, new_p1,new_p2, proj)
                        return newdeg
                        #save_degree(newdeg, newpA,newpB, new_p1, new_p2, nl1,nl2, method, frac,case, "after",proj)

                        #if method == "discrete":
                        #    old_deg1.append(deg)
                        #    new_deg1.append(newdeg)
                        #if method == "itplane":
                        #    old_deg2.append(deg)
                        #    new_deg2.append(newdeg)




def compute_curvature(dirpath, smooth, name, alpha, beta):
    case = dirpath[-5:]

    # Input filenames
    model_path = path.join(dirpath, name, "model.vtp")

    # Output names
    model_smoothed_path = path.join(dirpath, name, "model_smoothed.vtp")
    #model_new_surface = path.join(dirpath, name, "newmodel_alpha_%s_beta_%s_%s.vtp" % (alpha,beta, method))
    #model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
    #model_new_surface_tmp = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s_tmp.vtp" % (alpha, beta))
    dirpath2 = "modelstable"
    model_new_surface = path.join(dirpath2, "%s_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    new_centerlines_path = path.join(dirpath2, "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    new_centerlines_path = path.join(dirpath, name, "new_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta))
    new_centerlines_path = path.join(dirpath, name, "new_centerlines_alpha_0.2_beta_-0.1.vtp")
    sppath = path.join(dirpath, name, "centerline_clipped_part.vtp")

    # Centerliens
    centerline_complete_path = path.join(dirpath, name, "centerline_complete.vtp")
    centerline_clipped_path = path.join(dirpath, name,  "centerline_clipped.vtp")
    centerline_clipped_part_path = path.join(dirpath, name,  "centerline_clipped_part.vtp")
    centerline_clipped_part1_path = path.join(dirpath, name,  "centerline_clipped_end_part.vtp")
    centerline_new_path = path.join(dirpath, name, "centerline_interpolated.vtp")
    centerlines_ordered_path = path.join(dirpath, name, "centerline_ordered.vtp")

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
    
    # Compute centerline
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

    # Special cases including the opthalmic artery
    eye = False
    eye, clip_ID, centerlines_complete, eyeline = find_ophthalmic_artery(centerlines_complete, clipping_points)

    # Find longest line:
    lines = []
    n = centerlines_complete.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines_complete, i))

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
    line = ExtractSingleLine(centerlines_in_order, 0)
    print("Clipping centerlines.")
    patch_cl = CreateParentArteryPatches(centerlines_in_order, clip_points)
    
    # Get clipped curve 
    print("Getting clipped curve.")
    clipped_curve = get_clipped_curve(centerlines_in_order, clipping_points) 
    
    patch_start = ExtractSingleLine(patch_cl, 0)
    patch_ends = []
    for i in range(1,n+1):
        patch_ends.append(ExtractSingleLine(patch_cl, i))

    # Find ID of middle pooint:
    growDir_tmp = "horizont"
    middle_points, middleIds = get_spline_points(line, beta, growDir_tmp, dirpath, clip_points,voronoi)

    locator = get_locator(line)
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)
    dx_p1 = middle_points[0] - p1
    dx_p2 = middle_points[-1] - p2
    middle_points = middle_points[1:-1]

    curv_before =  [] 
    curv_after =  [] 

    # Compute new centerline using VMTK
    compute_vmtk = True
    if compute_vmtk:
        new_centerline_vmtk = centerlines_in_order#ReadPolyData(new_centerlines_path)#makeCenterline(model_new_surface, new_centerlines_path, smooth=False, resampling=False)
        #new_centerline_vmtk = ReadPolyData(new_centerlines_path)#makeCenterline(model_new_surface, new_centerlines_path, smooth=False, resampling=False)
        #new_cl =  makeCenterline(model_path, centerline_complete_path, smooth=False, resampling=False)
        new_cl = centerlines_in_order
        
        WritePolyData(new_cl, "oldcl.vtp")
        WritePolyData(new_centerline_vmtk, "newcl.vtp")

        lines = []
        n = new_centerline_vmtk.GetNumberOfCells()
        for i in range(n):
            lines.append(ExtractSingleLine(new_centerline_vmtk, i))

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

        new_line_vmtk = longest[0]
        curvs = []

        locator = get_locator(line)
        IDP0_1 = locator.FindClosestPoint(p1)
        IDPN_1 = locator.FindClosestPoint(p2)

        locator = get_locator(new_line_vmtk)
        IDP0_2 = locator.FindClosestPoint(p1)
        IDPN_2 = locator.FindClosestPoint(p2)

        plt.axvline(x=IDP0_1,ls='--')
        plt.axvline(x=IDPN_1,ls='--')

        # 1) VMTK - Factor variance
        factors = [0.2, 0.5,1.0,1.5,1.7,1.9]
        for factor in factors:
            factor = 0.7
            it = 100
            line_fac = CenterlineAttribiutes(new_line_vmtk, smooth=True, it=100, factor=factor)
            line_old = CenterlineAttribiutes(line, smooth=True, it=it, factor=factor)
            curv_fac = get_array("Curvature", line_fac)
            curv_old = get_array("Curvature", line_old)
            curv_fac= resample(curv_fac,curv_old.size)
            curv_fac= gauss(curv_fac, 5)
            curv_old= gauss(curv_old, 5)
            curvs.append(curv_fac)

            curv_before.append(max(curv_old[IDP0_1:IDPN_1])[0])
            curv_after.append(max(curv_fac[IDP0_2:IDPN_2])[0])
            embed()
            

        # 2) VMTK - Iteration variance
        its = [25,50,75,100,150,200]
        #its = [200]
        WritePolyData(new_line_vmtk,"newlinevmtk.vtp")

        for it in its:
            line_it = CenterlineAttribiutes(new_line_vmtk, smooth=True, it=it, factor=1.0)
            curv_it = get_array("Curvature", line_it)
            curv_it = resample(curv_it,1310 )
            curv_it= gauss(curv_it, 5)
            curvs.append(curv_it)

            line_old = CenterlineAttribiutes(line, smooth=True, it=it, factor=1.0)
            curv_old = get_array("Curvature", line_old)
            curv_old= gauss(curv_old, 5)

            curv_before.append(max(curv_old[IDP0_1:IDPN_1])[0])
            curv_after.append(max(curv_fac[IDP0_2:IDPN_2])[0])


        for i,c in enumerate(curvs):
            if i < len(factors):
#
                #plt.plot(c, label="VMTK, Factor=%.1f" % (factors[i]), linewidth=2.5)
                plt.plot(c, label="Factor=%.1f" % (factors[i]), linewidth=2.5)
                pass
            else:
                #plt.plot(c, label="Iterations=%i" % (its[i-len(factors)]), linewidth=2.5)
                #plt.plot(c, label="VMTK, Iterations=%i" % (its[i-len(factors)]), linewidth=2.5)
                pass

        ax = plt.gca()
        cl = True
        if cl:
            plt.xlabel(r"Abscissa", fontsize=20)
            plt.ylabel(r"Curvature",fontsize=20)
            plt.xticks(fontsize=15)
            #ax.xaxis.set_label_coords(0.85, -0.075)
            plt.yticks(fontsize=15)
            plt.legend(fancybox=True, shadow=True, fontsize=15, loc="upper center", bbox_to_anchor=(0.5,1.15), ncol=6)
        else:
            plt.xticks([])#fontsize=15)
            plt.yticks([])#fontsize=15)
        plt.show()

    # Compute new centerline by manual displacement
    print("Moving centerline manually")
    p_1 = patch_start
    p_2 = patch_ends[0]
    patch_cl_new1 = move_line_horo(p_1, ID1, ID2, dx_p1, clip=False, side="left")
    patch_cl_new2 = move_line_horo(p_2, ID1, ID2, dx_p1, clip=False, side="right")
    clipped_part_new = move_line_horo(clipped_curve, ID1, ID2, dx_p1, clip=True, eye=eye)
    new_centerline = merge_data([patch_cl_new1, clipped_part_new, patch_cl_new2])
    clipped_part_new = connect_line(clipped_part_new)
    growDir = "vertical"
    middle_points, middleIds, dx = get_spline_points(new_centerline, alpha, growDir, dirpath, clip_points,voronoi)
    clipped_part_new = move_dots_vertical(clipped_part_new, dx, eye=eye)
    new_centerline = merge_data([patch_cl_new1, clipped_part_new, patch_cl_new2])
    new_centerline = connect_line(new_centerline)
    

    locator = get_locator(line)
    IDP0_a = locator.FindClosestPoint(p1)
    IDPN_a = locator.FindClosestPoint(p2)
    
    locator = get_locator(new_centerline)
    IDP0_b = locator.FindClosestPoint(p1)
    IDPN_b = locator.FindClosestPoint(p2)

    # 3) Discrete derivatives
    #new_centerline_resamp = CenterlineResampling(new_centerline, 0.01, filename=None)
    disc = False
    if disc:
        neigh = 20
        #ns = [15,17,20,24,26,28,30,32]
        curvs = []
        ns = [10,15,20,25,30,35]
        ns = [20]
        #ax = plt.gca()
        #for neigh in ns:
        for i,neigh in enumerate(ns):
            #line_di_old, curv_di_old = discrete_geometry(dirpath, line, neigh=neigh)
            # KEEP THIS LINE FOR true
            line_di, curv_di = discrete_geometry(dirpath, new_centerline, neigh=neigh)

           #     print line_di.GetNumberOfPoints()
            filtercurv = gauss(curv_di,5)
            curvs.append(filtercurv)
            siph = ReadPolyData(sppath)
            ll = ExtractSingleLine(siph,0)
            WritePolyData(line, "helphelp.vtp")
            WritePolyData(ll, "helpSIPHONhelp.vtp")
            embed()

        for i,filtercurv in enumerate(curvs):
        #plt.plot(curv_di,label="no resamp")
        #plt.plot(curv_di_0,label="initial")
            plt.plot(filtercurv, label=r"$m$=%s" % ns[i], linewidth=2.5)
            #plt.plot(filtercurv, label=r"Discrete derivative, $m$=%s" % ns[i], linewidth=2.5)
            #maxcurv = max(filtercurv[IDP0_b+10:IDPN_b-10])
            #return maxcurv

            #line_old, curv_old = discrete_geometry(dirpath, line, neigh=neigh)
            #curv_old= gauss(curv_old, 5)

            #curv_before.append(max(curv_old[IDP0_a:IDPN_a]))
            #curv_after.append(max(curv_di[IDP0_b:IDPN_b]))


        cl = True
        if cl:
            plt.legend(fancybox=True, shadow=True, fontsize=15, loc="upper center", bbox_to_anchor=(0.5,1.15), ncol=len(ns))
            plt.xlabel(r"Abscissa", fontsize=20)
            plt.ylabel(r"Curvature",fontsize=20)
            plt.xticks(fontsize=15)
            #ax.xaxis.set_label_coords(0.85, -0.075)
            plt.yticks(fontsize=15)
        else:
            plt.xticks([])#fontsize=15)
            plt.yticks([])#fontsize=15)


        plt.axvline(x=IDP0_b,ls='--')
        plt.axvline(x=IDPN_b,ls='--')

        plt.show()

   # return maxcurv
    #print IDP0
    #print IDPN
    #plt.axvline(x=IDP0,ls='--')
    #plt.axvline(x=IDPN,ls='--')
    if False:
        plt.title(r"Case %s, $\alpha$=%s $\beta$=%s" % (dirpath[-5:], alpha,beta))
        plt.plot(curv_di, label="After")
        plt.legend()

        plt.axvline(x=IDP0,ls='--')
        plt.axvline(x=IDPN,ls='--')
        plt.show()
    #plt.show()
    #return maxcurv

    # 4) Knot free regression splines (Discrete diff init)
    discrete_init = False
    if discrete_init:
        clfile = writeCenterline(new_centerline)
        curv_di_cut = []
        for c in curv_di: # Ignore outer points
            if c > 0.001:
                curv_di_cut.append(c)
        order = 5
        curv_new = np.array(curv_di_cut)
        yhat = savitzky_golay(curv_new, 51, 5)
        start = [0 for i in range(order)]
        endp = get_curvilinear_coordinate(new_centerline)[-1]
        end = [endp for i in range(order)]
        length = get_curvilinear_coordinate(new_centerline)
        curv_newnew = resample(yhat, len(length))
        top = list(argrelextrema(curv_newnew, np.greater)[0])
        bot = list(argrelextrema(curv_newnew, np.less)[0])
        minmax = sorted(top + bot)
        middle = []
        for ID in minmax:
            middle.append(length[ID])
        
        initial_knot_position = start + middle + end
        curv_kf_discinit = spline_matlab(".", clfile, init_knots=20, order=float(order),plot=False, init_array=initial_knot_position)

    # 5) Knot free regression splines 
    knot_free = False
    if knot_free:
        clfile = writeCenterline(new_centerline)
        kfcurvs = []
        inits = [10,15,20,25,30,35,40]
        #inits = [50 + 2*i for i in range(20)]
        #inits = [20]
        ax = plt.gca()

        for I in inits:

            #curv_kf = spline_matlab(".", clfile, init_knots=I, order=5,plot=False, init_array=initial_knot_position)
            curv_kf = spline_matlab(".", clfile, init_knots=I, order=float(5),plot=False)
            curv_kf = gauss(curv_kf,5)
            kfcurvs.append(curv_kf)
            curv_after.append(max(curv_kf[IDP0_b:IDPN_b]))

        clfile = writeCenterline(line)
        for I in inits:
            curv_old = spline_matlab(".", clfile, init_knots=I, order=float(5),plot=False)
            curv_old= gauss(curv_old,5)
            curv_before.append(max(curv_old[IDP0_a:IDPN_a]))

        for i, c in enumerate(kfcurvs):
            plt.plot(c, label=r"Init knots=%i" % inits[i], linewidth=2.5)
            #plt.plot(c, label=r"Free-knot Spline, Init knots=%i" % inits[i], linewidth=2.5)

        cl = True
        if cl:
            plt.xlabel(r"Abscissa", fontsize=20)
            plt.ylabel(r"Curvature",fontsize=20)
            plt.xticks(fontsize=15)
            #ax.xaxis.set_label_coords(0.85, -0.075)
            plt.yticks(fontsize=15)
            plt.legend(fancybox=True, shadow=True, fontsize=14, loc="upper center", bbox_to_anchor=(0.5,1.15), ncol=len(inits))
        else:
            plt.xticks([])#fontsize=15)
            plt.yticks([])#fontsize=15)

        plt.show()


    # 6) Splines - 50 knots
    splines=False
    if splines:
        nknots = 50

        nk = [46,47,48,49,50]
        #nk = [55 + i for i in range(25)]
        #nk = [12 + 4*i for i in range(10)]
        nk = [10,15,20,25,30,35]
        #nk = [20]
        nkcurvs = []
        ax = plt.gca()
        for nknots in nk:
            line_sp, curv_sp = spline_line(new_centerline, get_curv=True, isline=True, nknots=nknots)
            curv_sp = gauss(curv_sp,5)
            #print max(curv_sp[IDP0:IDPN])
            nkcurvs.append(curv_sp)

            line_old, curv_old = spline_line(line, get_curv=True, isline=True, nknots=nknots)
            curv_old= gauss(curv_old,5)
            curv_before.append(max(curv_old[IDP0_a:IDPN_a]))
            curv_after.append(max(curv_sp[IDP0_b:IDPN_b]))

        for i, c in enumerate(nkcurvs):
            plt.plot(c, label=r"Knots=%i" % nk[i], linewidth=2.5)
          #  plt.plot(c, label=r"B-Spline, Knots=%i" % nk[i], linewidth=2.5)

        cl = True
        if cl:
            plt.xlabel(r"Abscissa", fontsize=20)
            plt.ylabel(r"Curvature",fontsize=20)
            plt.xticks(fontsize=15)
            #ax.xaxis.set_label_coords(0.85, -0.075)
            plt.yticks(fontsize=15)
            plt.legend(fancybox=True, shadow=True, fontsize=14, loc="upper center", bbox_to_anchor=(0.5,1.15), ncol=len(nk))
        else:
            plt.xticks([])#fontsize=15)
            plt.yticks([])#fontsize=15)

        plt.show()


    print len(curv_di)
    curv_fac = resample(curv_fac, len(curv_di))
    curv_it = resample(curv_it, len(curv_di))
    curv_fac = gauss(curv_fac, 5)
    curv_it = gauss(curv_it, 5)
    #plt.plot(curv_fac, '-',color="brown", label=r"VMTK, Factor=1.5" , linewidth=2.5)
    #plt.plot(curv_it,'--', color="brown", label=r"VMTK, Iterations=200", linewidth=2.5)
    #plt.legend(fancybox=True, shadow=True, fontsize=14, loc="upper center", bbox_to_anchor=(0.5,1.15), ncol=5)

    #plt.show()
    #plt.legend(fancybox=True, shadow=True, fontsize=15)
    #plt.show()

    # Save curvatures 
    print curv_before
    print curv_after
    methods = ["fac", "it", "discrete", "freeknot", "splines"]
    params = [factors[0], its[0], ns[0], inits[0], nk[0]]
    times = ["before", "after"]
    curv_array = [curv_before, curv_after]
    for j,curvs in enumerate(curv_array):
        time = times[j]
        for i, curv in enumerate(curvs):
            method = methods[i]
            param = params[i]
  #          save_curv(curv, method, param ,case, time )

    """
    locator = get_locator(new_centerline_resamp)
    IDP0 = locator.FindClosestPoint(p1)
    IDPN = locator.FindClosestPoint(p2)
    clfile = writeCenterline(new_centerline_resamp)


    curv_kf3 = spline_matlab(".", clfile, init_knots=20, order= 5.,plot=False)
    curv_kf4 = spline_matlab(".", clfile, init_knots=20, order= 6.,plot=False)
    
    ### Compare initial knots ###
    m_c = []
    inits=[15,16,17,18,19]
    for init in inits: 
        curv_m = spline_matlab(".", clfile, init_knots=init, order= 5.,plot=False)
        m_c.append(curv_m)
    for i, c in enumerate(m_c):
        plt.plot(c, label="Init_knots=%s" %inits[i] ,linewidth=2)
    plt.xlabel("Centerline length, resampled")
    plt.ylabel("Curvature")
    plt.axvline(x=IDP0,ls='--')
    plt.axvline(x=IDPN,ls='--')
    plt.legend()
    plt.show()

    clfile = writeCenterline(new_centerline)
    inits=[15,16,17,18,19]
    m_c = []
    locator = get_locator(new_centerline)
    IDP0 = locator.FindClosestPoint(p1)
    IDPN = locator.FindClosestPoint(p2)
    for init in inits: 
        curv_m = spline_matlab(".", clfile, init_knots=init, order= 5.,plot=False)
        m_c.append(curv_m)
    for i, c in enumerate(m_c):
        plt.plot(c, label="Init_knots=%s" %inits[i] ,linewidth=2)
    plt.xlabel("Centerline length")
    plt.ylabel("Curvature")
    plt.axvline(x=IDP0,ls='--')
    plt.axvline(x=IDPN,ls='--')
    plt.legend()
    plt.show()

    # Discard ends of line for comparisson
    ID0 = 50
    IDN = len(curv_sp) - 50
    
    P0 = line_sp.GetPoint(ID0)
    PN = line_sp.GetPoint(IDN)

    locator = get_locator(new_line_vmtk)
    IDP0 = locator.FindClosestPoint(P0)
    IDPN = locator.FindClosestPoint(PN)


    # Truncate curvature arrays
    curv_fac = curv_fac[IDP0 : IDPN]
    curv_it  = curv_it[IDP0 : IDPN]
    curv_di  = curv_di[ID0:IDN]
    curv_sp  = curv_sp[ID0:IDN]



    locator = get_locator(new_centerline_resamp)
    IDA = locator.FindClosestPoint(P0)
    IDB = locator.FindClosestPoint(PN)

    curv_kf  = curv_kf[IDA:IDB]
    #curv_kf2  = curv_kf2[IDA:IDB]




    #curv_kf3  = curv_kf3[IDA:IDB]
    #curv_kf4  = curv_kf4[IDA:IDB]

    # Make all lines equal length
    curv_fac = resample(curv_fac, len(curv_di))
    curv_it  = resample(curv_it , len(curv_di))

    curv_kf  = resample(curv_kf , len(curv_di))
    #curv_kf2  = resample(curv_kf2 , len(curv_di))

    #curv_kf3  = resample(curv_kf3, len(curv_kf))
    #curv_kf4 = resample(curv_kf4, len(curv_kf))

    # Find clipping points
    locator = get_locator(line_sp)
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)

    # Plot
    curvs = [curv_kf, curv_fac, curv_it, curv_sp, curv_di]
    labels = ["Knot free, initial knots = %s, order = %s", "Factor = %.1f", "Iterations = %s", "Splined, knots = %s", "Discrete, neighbours = %s"]
    values = [20, 1.0, 100, 50, 10]

    #plt.plot(curv_kf2, 'r--', label=labels[0] % (values[0],6) , linewidth=2)
    #plt.plot(curv_kf3, 'r:', label=labels[0] % (values[0],5) + " , resampled" , linewidth=2)
    #plt.plot(curv_kf4, 'r-.', label=labels[0] % (values[0],6) + " ,resampled", linewidth=2)

    plt.plot(curvs[0], 'r-', label=labels[0] % (values[0],5) , linewidth=2)
    plt.plot(curvs[1], 'y--', label=labels[1] % values[1] ,linewidth=2)
    plt.plot(curvs[2], 'y-', label=labels[2] % values[2], linewidth=2)
    plt.plot(curvs[3], 'g-', label=labels[3] % values[3], linewidth=2)
    plt.plot(curvs[4], 'm-', label=labels[4] % values[4], linewidth=2)

    plt.legend()
    plt.xlabel("Centerline length")
    plt.ylabel("Curvature")
    plt.axvline(x=ID1,ls='--')
    plt.axvline(x=ID2,ls='--')
    plt.show()
    """
    return 0#maxcurv


            #newp2_2, np2_id = find_closest_point(T2[-1],newmaxID, moved_siphon.GetNumberOfPoints(), moved_siphon, show=False)
def find_closest_point(dx,start,stop, P0, line, show=True):

    a = dx[0]
    b =  dx[1]
    c =  dx[2]
    n = np.array([a,b,c])
    n = n / la.norm(n)

    # Define plane
    xmin = 42#min(Z, key=operator.itemgetter(1))[0] - 4 
    xmax = 52#max(Z, key=operator.itemgetter(1))[0] + 4
    ymin = 20#min(Z, key=operator.itemgetter(1))[1] - 4 
    ymax = 29#max(Z, key=operator.itemgetter(1))[1] + 4
    xx,yy = np.meshgrid(np.linspace(xmin,xmax,150),np.linspace(ymin,ymax,150))
    d = a*P0[0] + b*P0[1] + c*P0[2]  
    zz = (d - a*xx - b*yy ) / float(c)
    points = []
    for i in range(start,stop):
        p = line.GetPoint(i)
        points.append( np.array(p) )
    dist_list = []
    for i,pcl in enumerate(points):
        v = pcl - np.array(P0)
        dist = abs(v.dot(n))
        dist_list.append(dist)

    
    minID = dist_list.index(min(dist_list)) + start
    minP = points[minID - start]
    if show:
        plotPlane(xx,yy,zz, points, minP)
    return minP, minID 

def find_furthest_points(dx,p1 ,p2 , line, show=True):
    P0 = line.GetPoint(0)
    a = dx[0]
    b =  dx[1]
    c =  dx[2]
    n = np.array([a,b,c])
    n = n / la.norm(n)

    # Define plane
    xmin = 42#min(Z, key=operator.itemgetter(1))[0] - 4 
    xmax = 52#max(Z, key=operator.itemgetter(1))[0] + 4
    ymin = 20#min(Z, key=operator.itemgetter(1))[1] - 4 
    ymax = 29#max(Z, key=operator.itemgetter(1))[1] + 4
    xx,yy = np.meshgrid(np.linspace(xmin,xmax,150),np.linspace(ymin,ymax,150))
    d = a*P0[0] + b*P0[1] + c*P0[2]  
    zz = (d - a*xx - b*yy ) / float(c)
    points = []
    for i in range(line.GetNumberOfPoints()):
        p = line.GetPoint(i)
        points.append(np.array(p))
    dist_list = []
    for i,pcl in enumerate(points):
        v = pcl - np.array(P0)
        dist = abs(v.dot(n))
        dist_list.append(dist)

    
    maxID = dist_list.index(max(dist_list))
    maxP = points[maxID]
    if show:
        plotPlane(xx,yy,zz, points, maxP)
    return maxP, maxID 


def plotPlane(xx,yy,zz,P, maxP):
    # Plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(xx,yy,zz, alpha=0.2)

    for x_,y_,z_ in P:
        ax.scatter(x_,y_,z_, color="blue")

    
    ax.scatter(maxP[0],maxP[1],maxP[2] ,color="red",  s=100)

    plt.show()


def spline_matlab(path, filename, init_knots, order, plot=True):
    """
    Perform Knot-free regresion spline on input centerline
    extracted from input filename
    """

    print("Computing knot free regression spline")
    #init_array_matlab = matlab.double( [i for i in init_array] )
    curv_m = mlab.CenterlineCharacterization(path, filename, init_knots, order, nargout=1)
    #curv_m = mlab.CenterlineCharacterization(path, filename, init_knots, order, init_array_matlab, nargout=1)
    n = len(curv_m)
    curv_p = np.zeros(n)

    for i in range(n):
        curv_p[i] = curv_m[i][0]

    if plot: 
        plt.plot(curv_p, label="Knot free regression spline")
        plt.legend()
        plt.show()

    return curv_p

def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    import numpy as np
    from math import factorial
       
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if  window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
       # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
       # pad the signal at the extremes with
       # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def odr_line(ID1, ID2, line, curvature, limit, boxplot=False):
    if limit == "cumulative":
        boxplot = False
    else:
        boxplot = True

    lim = len(curvature)-1 
    if not boxplot:
        max_cum = 10
        ID1_up = ID1+ 1
        ID1_down = ID1- 1
        ID2_up = ID2 - 1
        ID2_down = ID2 + 1
        while sum(curvature[ID1:ID1_up+1]) < max_cum and ID1_up < lim:
            ID1_up += 1
        while sum(curvature[ID1_down:ID1+1]) < max_cum and ID1_down > 0 :
            ID1_down -= 1
        while sum(curvature[ID2_up:ID2+1]) < max_cum  and ID2_up > 0:
            ID2_up -= 1
        while sum(curvature[ID2:ID2_down+1]) < max_cum and ID2_down < lim:
            ID2_down += 1
    else:
        SD = 0.045
        ID1_up = ID1 + 5
        ID1_down = ID1- 5
        ID2_up = ID2 - 5
        ID2_down = ID2 + 5
        mean1 = sum(curvature[ID1_down:ID1_up+1]) / 11.
        mean2 = sum(curvature[ID2_up:ID2_down+1]) / 11.
        SD1 = np.sqrt( sum( (curvature[ID1_down:ID1_up+1] - mean1)**2) / 10 )
        SD2 = np.sqrt( sum( (curvature[ID2_up:ID1_down+1] - mean2)**2) / 10 )
        tol1 = mean1 + SD1*1.96
        tol2 = mean2 + SD2*1.96
        while curvature[ID1_up] < tol1 and ID1_up < lim:
  #          print curvature[ID1_up]
            ID1_up +=1
        while curvature[ID1_down] < tol1 and ID1_down > 0:
            ID1_down -=1
        while curvature[ID2_up] < tol2 and ID2_up > 0:
            ID2_up -=1
        while curvature[ID2_down] < tol2 and ID2_down < lim:
            ID2_down +=1

    p1s = []
    for i in range(ID1_down, ID1_up+1):
        p1s.append(line.GetPoint(i))

    p2s = []
    for i in range(ID2_up, ID2_down+1):
        p2s.append(line.GetPoint(i))

    # Arrange points in matrix
    X1 = np.array([list(p) for p in p1s])
    X2 = np.array([list(p) for p in p2s])
    
    # Find mean of points
    avg1 = np.array([np.mean(X1[:,0]), np.mean(X1[:,1]), np.mean(X1[:,2])])        
    avg2 = np.array([np.mean(X2[:,0]), np.mean(X2[:,1]), np.mean(X2[:,2])])        
    # Subtract the mean from all points
    dX1 = X1 - np.array([avg1 for i in range(len(X1))])
    dX2 = X2 - np.array([avg2 for i in range(len(X2))])

    # Find SVD
    U,S,V1 = la.svd(dX1)
    U,S,V2 = la.svd(dX2)
    # Find direction vector
    d1 = V1[0]
    d2 = V2[0]
    p01 = avg1
    p02 = avg2
    # Parametric equation P = p0 + t*d
    # Make lines with curv
    # Create edges between new_centerline points
    curv_lines_split = []
    points = [p1s,p2s]
    for k, p in enumerate(points):
        pts = vtk.vtkPoints()
        for i in range(len(p)):
            pts.InsertNextPoint(p[i])

        lines = vtk.vtkCellArray()
        for i in range(len(p)-2):
            newline = vtk.vtkLine()
            newline.GetPointIds().SetId(0, i)
            newline.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(newline)

        line_ = vtk.vtkPolyData()
        line_.SetPoints(pts)
        line_.SetLines(lines)

        M = line_.GetNumberOfPoints()
        curvArray = get_vtk_array("Curvature", 1 ,M)
        if k == 0:
            for i in range( ID1_up+1- ID1_down):
                curvArray.SetTuple(i, [curvature[ID1_down + i]])
        else:
            for i in range( ID2_down+1 - ID2_up):
                curvArray.SetTuple(i, [curvature[ID2_up + i]])

        line_.GetPointData().AddArray(curvArray)

        curv_lines_split.append(line_)

    curvlines = merge_data(curv_lines_split)
    return p01,d1, p02, d2, curvlines


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


def discrete_geometry(case, line, neigh=10):
    len_line = line.GetNumberOfPoints()
    N = line.GetNumberOfPoints() 

    # Compute cumulative chord length
    t = np.zeros(N)
    p = []
    for i in range(N):
        p.append(np.array(list(line.GetPoint(i))))
        p[i] = np.array(p[i])

    norms = [la.norm(p[j] - p[j-1]) for j in range(1,N)]
    s = sum(norms)
    for i in range(1, N):
        s1 = sum(norms[:i+1])
        t[i] = s1 / s

    # Radius of sliding neighbourhood
    m = neigh

    dxdt = np.zeros(N)
    dydt = np.zeros(N)
    dzdt = np.zeros(N)

    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)

    for i in range(N):
        x[i] = p[i][0]
        y[i] = p[i][1]
        z[i] = p[i][2]
 

    for i in range(0, m):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(0, 2*m+1)])
        dxdt[i] = sum([ (t[j] - t[i])*(x[j]-x[i]) for j in range(0, 2*m+1)]) / t_sum
        dydt[i] = sum([ (t[j] - t[i])*(y[j]-y[i]) for j in range(0, 2*m+1)]) / t_sum
        dzdt[i] = sum([ (t[j] - t[i])*(z[j]-z[i]) for j in range(0, 2*m+1)]) / t_sum

    for i in range(m, N-m):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(i-m, i+m+1)])
        dxdt[i] = sum([ (t[j] - t[i])*(x[j]-x[i]) for j in range(i-m, i+m+1)]) / t_sum
        dydt[i] = sum([ (t[j] - t[i])*(y[j]-y[i]) for j in range(i-m, i+m+1)]) / t_sum
        dzdt[i] = sum([ (t[j] - t[i])*(z[j]-z[i]) for j in range(i-m, i+m+1)]) / t_sum

    for i in range(N-m, N):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(N-2*m, N)])
        dxdt[i] = sum([ (t[j] - t[i])*(x[j]-x[i]) for j in range(N-2*m-1, N)]) / t_sum
        dydt[i] = sum([ (t[j] - t[i])*(y[j]-y[i]) for j in range(N-2*m-1, N)]) / t_sum
        dzdt[i] = sum([ (t[j] - t[i])*(z[j]-z[i]) for j in range(N-2*m-1, N)]) / t_sum

    dgammadt = []
    dgammadt_norm = np.zeros(N)
    for i in range(N):
        dgammadt.append(np.array([dxdt[i], dydt[i], dzdt[i]]))
        dgammadt_norm[i] = la.norm(dgammadt[i])
    
    tg = []
    for i in range(N):
        tg.append(dgammadt[i] / dgammadt_norm[i])

    t1 = np.zeros(N)
    t2 = np.zeros(N)
    t3 = np.zeros(N)

    for i in range(N):
        t1[i] = tg[i][0]
        t2[i] = tg[i][1]
        t3[i] = tg[i][2]
    
    dt1dt = np.zeros(N)
    dt2dt = np.zeros(N)
    dt3dt = np.zeros(N)
    #m #= m+1


    for i in range(0, m):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(0, 2*m+1)])
        dt1dt[i] = sum([ (t[j] - t[i])*(t1[j]-t1[i]) for j in range(0, 2*m+1)]) / t_sum
        dt2dt[i] = sum([ (t[j] - t[i])*(t2[j]-t2[i]) for j in range(0, 2*m+1)]) / t_sum
        dt3dt[i] = sum([ (t[j] - t[i])*(t3[j]-t3[i]) for j in range(0, 2*m+1)]) / t_sum

    for i in range(m, N-m):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(i-m, i+m+1)])
        dt1dt[i] = sum([ (t[j] - t[i])*(t1[j]-t1[i]) for j in range(i-m, i+m+1)]) / t_sum
        dt2dt[i] = sum([ (t[j] - t[i])*(t2[j]-t2[i]) for j in range(i-m, i+m+1)]) / t_sum
        dt3dt[i] = sum([ (t[j] - t[i])*(t3[j]-t3[i]) for j in range(i-m, i+m+1)]) / t_sum
    
    for i in range(N-m, N):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(N-2*m, N)])
        dt1dt[i] = sum([ (t[j] - t[i])*(t1[j]-t1[i]) for j in range(N-2*m-1, N)]) / t_sum
        dt2dt[i] = sum([ (t[j] - t[i])*(t2[j]-t2[i]) for j in range(N-2*m-1, N)]) / t_sum
        dt3dt[i] = sum([ (t[j] - t[i])*(t3[j]-t3[i]) for j in range(N-2*m-1, N)]) / t_sum

    dtgdt = []
    dtgdt_norm = np.zeros(N)
    for i in range(N):
        dtgdt.append(np.array([dt1dt[i], dt2dt[i], dt3dt[i]]))
        dtgdt_norm[i] = la.norm(dtgdt[i])

    curv = np.zeros(N)
    for i in range(N):
        curv[i] = dtgdt_norm[i] / dgammadt_norm[i]

    curv = resample(curv, len_line)

    return line, curv

def spline_line(line, get_curv=False, isline=False, nknots=50):
    
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

def cutLine(line, endID):
    # Create edges between new_centerline points
    pts = vtk.vtkPoints()
    for i in range(endID):
        pts.InsertNextPoint(line.GetPoint(i))

    lines = vtk.vtkCellArray()
    for i in range(endID-2):
        newline = vtk.vtkLine()
        newline.GetPointIds().SetId(0, i)
        newline.GetPointIds().SetId(1, i + 1)
        lines.InsertNextCell(newline)

    line = vtk.vtkPolyData()
    line.SetPoints(pts)
    line.SetLines(lines)
    return line

def find_ophthalmic_artery(centerlines, clip_pts):

    # Extract lines:
    lines = []
    n = centerlines.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines, i))

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

def get_spline_points(line, param,growDir,case, clip_points, voronoi):
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
    
    # Select n uniformly spaced points 
    n = 7
    points = []
    ids = np.zeros(n)
    dx = 1 / (n + 1.)
    for i in range(1,n+1):
        ID = int(ID1 + (ID2 - ID1) * i * dx)
        ids[i-1] = ID
        p = line.GetPoints().GetPoint(ID)
        points.append(np.array([p[0], p[1], p[2]]))

    for i in range(len(P)):
        P[i] = np.array([P[i][0], P[i][1], P[i][2]])

    xx,yy,zz,n = bestPlane(points,P)

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

def writeCenterline(centerline):
    P = []
    for i in range(centerline.GetNumberOfPoints()):
        P.append(np.array(centerline.GetPoint(i)))

    clfile = "centerlines/clpoints_discrete.txt"
    with open(clfile,'wb') as f:
        for p in P:
            f.write("%s %s %s\n" % (p[0],p[1],p[2]))
    return clfile

def writeParticles(points, filename):
    if path.exists(filename):
        remove(filename)

    with open(filename, "a") as text:
        for P in points:
            text.write( "%.5f %.5f %.5f\n" % (P[0], P[1], P[2])  )

def writePoints(p1,p2,p3,p4, filename):
    if path.exists(filename):
        remove(filename)
        
    with open(filename, "a") as text:
        text.write( "%.5f %.5f %.5f\n" % (p1[0], p1[1], p1[2])  )
        text.write( "%.5f %.5f %.5f\n" % (p2[0], p2[1], p2[2])  )
        text.write( "%.5f %.5f %.5f\n" % (p3[0], p3[1], p3[2])  )
        text.write( "%.5f %.5f %.5f\n" % (p4[0], p4[1], p4[2])  )


def find_angle(pA, pB, p1,p2, proj):
    if not proj:
        P1A = np.array([pA[0] - p1[0], pA[1] - p1[1], pA[2] - p1[2]])
        P2B = np.array([pB[0] - p2[0], pB[1] - p2[1], pB[2] - p2[2]])
    else:
        P1A = np.array([0, pA[1] - p1[1], pA[2] - p1[2]])
        P2B = np.array([0, pB[1] - p2[1], pB[2] - p2[2]])
    costheta = (P1A.dot(P2B)) / (la.norm(P1A)*la.norm(P2B))
    angle = np.arccos(costheta)
    deg = (angle * 180 / np.pi)
    return deg, P1A, P2B


def find_angle_odr(d1,d2,proj):
    if d1.dot(d2) > 0:
        d1 = -d1

    if proj:
        d1[0] = 0
        d2[0] = 0

    costheta = (d1.dot(-d2)) / (la.norm(d1)*la.norm(-d2))
    angle = np.arccos(costheta)
    deg = (angle * 180 / np.pi)
    return deg, d1, d2

def print_curv_info():
    methods = ["fac", "it", "discrete", "freeknot", "splines"]
    params = [1.5, 200, 20, 20, 20]
    times = ["before", "after"]

    files = sorted(listdir("curvature"))
    before = files[1::2]
    after = files[::2]
    k=0
    for old, new in zip(before,after):
        # Extract angle values
        old_deg = open("curvature/" + old, "r")
        new_deg = open("curvature/" + new, "r")
        old_curv = []
        new_curv = []
        info = old.split("_")
        method = info[0]
        param = info[1]

        for i in range(10):
            oldinfo = old_deg.readline().split()
            newinfo = new_deg.readline().split()
            old_curv.append(float(oldinfo[-1]))
            new_curv.append(float(newinfo[-1]))

        # Print table values
        p0 = np.array(old_curv)
        p1 = np.array(new_curv)
        dS = p1- p0
        mean0 = sum(p0) / len(p0) 
        mean1 = sum(p1) / len(p1) 
        sd0 = np.sqrt(sum([ (xi - mean0)**2 for xi in p0]) / ( len(p0)-1))
        sd1 = np.sqrt(sum([ (xi - mean1)**2 for xi in p1]) / ( len(p1)-1))
        md = sum(dS) / len(dS)
        print "%s - %s & %.3f$\pm$%.3f & %.3f$\pm$%.3f & %.3f\\\\" % (method, param, mean0, sd0, mean1, sd1, md)

    return 

def print_angle_info():
    # Compute Sd, MEAN and statistics
    #old_deg = np.array(old_deg)
    #new_deg = np.array(new_deg)

    limits = ["cumulative","sd"] 

    methods = ["odrline", "MISR", "frac", "plane", "itplane", "itplane_clip", "maxcurv", "smooth", "discrete", "maxdist"] 
    files = sorted(listdir("angles"))
    before = files[1::2]
    after = files[::2]
    k=0
    for old, new in zip(before,after):
        # Extract angle values
        old_deg = open("angles/" + old, "r")
        new_deg = open("angles/" + new, "r")
        old_2D = []
        new_2D = []
        old_3D = []
        new_3D = []
        info = old.split("_")
        if len(info) < 4:
            method = info[0]
            param = info[1]
        else:
            method = info[0] + info[1]
            param = info[2]

        for i in range(20):
            oldinfo = old_deg.readline().split()
            newinfo = new_deg.readline().split()
            if i < 10:
                old_2D.append(float(oldinfo[-1]))
                new_2D.append(float(newinfo[-1]))
            else:
                old_3D.append(float(oldinfo[-1]))
                new_3D.append(float(newinfo[-1]))

        # Print table values
        tod = False
        threed = True
        if tod:
            if isnan(new_2D[-1]):
                p0 = np.array(old_2D[:-1])
                p1 = np.array(new_2D[:-1])
            else:
                p0 = np.array(old_2D)
                p1 = np.array(new_2D)
            dS = p1- p0
            mean0 = sum(p0) / len(p0) 
            mean1 = sum(p1) / len(p1) 
            sd0 = np.sqrt(sum([ (xi - mean0)**2 for xi in p0]) / ( len(p0)-1))
            sd1 = np.sqrt(sum([ (xi - mean1)**2 for xi in p1]) / ( len(p1)-1))
            md = sum(dS) / len(dS)
            print "%s - %s & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.2f\\\\" % (method, param, mean0, sd0, mean1, sd1, md)
            
        if threed:
            if isnan(new_3D[-1]):
                p0 = np.array(old_3D[:-1])
                p1 = np.array(new_3D[:-1])
            else:
                p0 = np.array(old_3D)
                p1 = np.array(new_3D)
            dS = p1- p0
            mean0 = sum(p0) / len(p0) 
            mean1 = sum(p1) / len(p1) 
            sd0 = np.sqrt(sum([ (xi - mean0)**2 for xi in p0]) / ( len(p0)-1))
            sd1 = np.sqrt(sum([ (xi - mean1)**2 for xi in p1]) / ( len(p1)-1))
            md = sum(dS) / len(dS)
            print "%s - %s& %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.2f\\\\" % (method, param, mean0, sd0, mean1, sd1, md)


    methods = [ "MISR"]
    multiples = [1,1.5,2,2.5]

    methods = [ "frac"]
    n = [3,5,5,7,7]
    l = [1,2,1,1,2]
    r = [2,3,4,6,5]
    #l[i]divn_values[i]

    methods = ["plane", "itplane", "itplane_clip", "maxcurv", "smooth", "discrete", "maxdist"] 

    l= [1 ,3,4]
    n= [2 ,5,5]



def save_degree(deg, pA,pB, p1, p2, l1,l2, method, param ,case, time, proj):
    filePoints = "angledata/%s_points_%s_%s_%s.particles" % (case, method, param, time)
    writePoints(pA,pB,p1,p2, filePoints)

    fileLine = "angledata/%s_line_%s_%s_%s.particles" % (case, method, param, time)
    pl = [np.array(p1) - 2*l1 + 0.01*i*l1 for i in range(1250)]
    pr = [np.array(p2) - 2*l2 + 0.01*i*l2 for i in range(1200)]
    points = pl+pr
    writeParticles(points,fileLine)

    fileDeg = "angles/%s_%s_%s.particles" % (method, param, time)
    with open(fileDeg, "a") as text:
        if proj:
            text.write( "2D %s %s\n" %  (case,  deg) )
        else:
            text.write( "3D %s %s\n" %  (case,  deg ))
            
def save_curv(curv, method, param ,case, time):
    fileCurv = "curvature/%s_%s_%s.particles" % (method, param, time)
    with open(fileCurv, "a") as text:
        text.write( "%s %s\n" %  (case,  curv) )
            
    




def writeToFile(values, case, param):
    mat = np.matrix(values)
    #with open('init_%s.txt' % (param),'wb') as f:
    with open('new_%s_%s.txt' % (param,case),'wb') as f:
        for line in mat:
            np.savetxt(f, line, fmt='%.3f')
        #for line in values:
            #np.savetxt(f, line, fmt="%.5f")
        #    f.write("%.5f\n" % line)

        #f.close()

def get_alphabeta(method):
    folder = "plots_%s"  % method
    filepath = path.join(folder, "optimal_alphabeta.txt")
    a= []
    b= []
    myfile = open(filepath,"r")
    for i in range(20):
        info = myfile.readline().split()
        alpha = float(info[-2].split("=")[1])
        beta = float(info[-1].split("=")[1])
        a.append(alpha)
        b.append(beta)
    return a,b




def main():
    smooth, basedir, case, kappa, theta, alphas, betas = read_command_line()
    #basedir = "."
    folders = sorted([folder for folder in listdir(basedir) if folder[:2] in ["P0"]])
    name = "surface"
    #folders = [folders[0], folders[3], folders[5]]

    n = 30
    ab_bound = np.loadtxt("alphabeta_bound.txt")
    # Curavture
    test=False
    printinfo = False
    tables = False
    method = "angle"
    #method = "curv"
    a,b = get_alphabeta(method)
    
    if tables:
        k= 0
        proj = False
        plus = []
        minus = []
        inits = []
        for folder in folders:
            for i in range(2):
                alpha = a[k]
                beta = b[k]
                print alpha,beta
                case = path.join(basedir,folder)
                maxcurv = compute_curvature(case ,smooth, name,  alpha,beta )
                #angle = compute_angle(case ,smooth, name,alpha, beta, proj)
                k +=1
                if i % 2 == 0:
                    minus.append(maxcurv)
                else:
                    plus.append(maxcurv)
        for folder in folders:
            case = path.join(basedir,folder)
            alpha = beta = 0.00000001
            maxcurv = compute_curvature(case ,smooth, name,  alpha,beta )
            inits.append(maxcurv)
        print plus
        print minus
        print inits

    if test:

        alpha_ = [0.0001]
        beta_ = [0.5]
        proj = False# proj = False = 3D angle

        #beta_ = [0.001]
        #beta_ = [0.001]
        #alpha_ = [0.5]
        #alpha_ = [0.5]
        for alpha,beta in zip(alpha_,beta_):
            for folder in folders: 
                case = path.join(basedir,folder)
                #maxcurv = compute_curvature(case ,smooth, name,  alpha,beta )
                angle = compute_angle(case ,smooth, name,alpha, beta, proj)

        """
        print old_deg1
        print new_deg1
        print old_deg2
        print new_deg2
        print old_deg3
        print new_deg3
        """
        case = sorted(listdir("cases"))
        news = [new_deg1, new_deg2, new_deg3]
        olds = [old_deg1, old_deg2, old_deg3]
        for k,new_deg in enumerate(news):
            for i,p in enumerate(new_deg):
                dA = p - olds[k][i]
                print"%s & %.3f & %.3f & %.3f\\\\" % (case[i], olds[k][i], p,dA)
                print ""

    if printinfo:
        print_angle_info()
     ##   print_curv_info()

    if kappa:
        init_curv = []
        it = 1
        k = 0
        for folder in folders:
            max_curv_values = np.zeros((n,n))
            amin,amax,bmin,bmax = ab_bound[k][0], ab_bound[k][1], ab_bound[k][2], ab_bound[k][3]

            alpha_ = np.linspace(amin,amax, n)
            beta_ = np.linspace(bmin,bmax,n)
            print("Working on case " + folder)
            case = path.join(basedir,folder)
            for i,alpha in enumerate(alpha_):
                for j,beta in enumerate(beta_):
                    print("Iteration %i of %i" % (it, n*n*10))
                    maxcurv = compute_curvature(case ,smooth, name,  alpha,beta )
                    max_curv_values[i,j] = maxcurv
                    #init_curv.append(maxcurv)
                    it +=1

            writeToFile(max_curv_values, folder, "curvature")
            k += 1
        #writeToFile(init_curv, folder, "curvature")
    # Angle
    if theta:
        it = 1#1#1#1#1#1#1#1#1#1#1#3000#1
        k = 0
        for folder in folders:
            angle_values = np.zeros((n,n))
            amin,amax,bmin,bmax = ab_bound[k][0], ab_bound[k][1], ab_bound[k][2], ab_bound[k][3]
            alpha_ = np.linspace(amin,amax, n)
            beta_ = np.linspace(bmin,bmax,n)
            print("Working on case "+ folder)
            case = path.join(basedir,folder)

            for i,alpha in enumerate(alpha_):
                for j,beta in enumerate(beta_):
                    print("Iteration %i of %i" % (it, n*n*10))
                    angle = compute_angle(case ,smooth, name,alpha, beta)
                    angle_values[i,j] = angle
                    it +=1

       #     writeToFile(angle_values, folder, "angle")
            k+= 1

 #       print_angle_info(old_deg,new_deg)

if  __name__ == "__main__":
    print("Starting matlab..")
    #mlab = matlab.engine.start_matlab()
    main()


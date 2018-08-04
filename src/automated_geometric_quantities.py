from argparse import ArgumentParser
from os import path, remove, listdir
from scipy.signal import argrelextrema, resample
from scipy.ndimage.filters import gaussian_filter as gauss

import matlab.engine
import operator
import sys
import numpy.linalg as la

# Local import
from common import *
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from moveandmanipulatetools import *

def read_command_line():
    """
    Read arguments from commandline
    """
    parser = ArgumentParser()

    parser.add_argument('-d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('-c', '--case', type=str, default=None, help="Choose case")
    parser.add_argument('-k', '--curvature', type=bool, default=False,
                        help="Compute curvature variation", metavar="curvature")
    parser.add_argument('-t', '--angle', type=bool, default=False,
                        help="Compute angle variation", metavar="angle")
    parser.add_argument("-a", "--alpha", type=float, default=None,
                        help="Compression factor in vertical direction, ranging from -1.0 to 1.0")
    parser.add_argument("-b", "--beta", type=float, default=None,
                        help="Compression factor in horizontal direction, ranging from -1.0 to 1.0")
    parser.add_argument('-mc', '--method_curv', type=str, default="disc",
            help="Method for computing curv. Available methods: disc | knotfree | vmtkfactor | vmtkit | spline")
    parser.add_argument('-ma', '--method_angle', type=str, default="plane",
            help="Method for computing siphon angle. Available methods: plane | itplane | itplane_clip | maxcurv | smooth | discrete | frac | odrline | MISR ")

    args = parser.parse_args()

    return args.dir_path, args.case, args.curvature, args.angle, args.alpha, args.beta, args.method_curv, args.method_angle


def compute_angle(dirpath, point_path, name, alpha,beta, method, proj=False):
    """
    Primary collection of methods for computing the angle of a vessel bend.
    Three main methods are currently implemented:
    1) ODR methods: odrline
    2) Tracing point methods: maxcurv, smooth, discrete, frac, MISR
    3) Relative tracing point methods: plane, itplane, itplane_clip

    Args:
        dirpath (str): Path to case folder.
        point_path (str): Filename of clipping point file.
        name (str): Directory where surface models are located.
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        method (str): Method used to compute angle.
        proj (bool): True / False for computing 2D / 3D angle.

    Returns:
        newdeg (float): New angle of a vessel bend from a manipulated centerline.
    """
    # Input filenames
    model_path = path.join(dirpath, name, "model.vtp")

    # Centerliens
    centerline_complete_path = path.join(dirpath, name, "centerline_complete.vtp")

    # Find endID from landmarking
    siphon_path = path.join(dirpath, "surface", "carotid_siphon.vtp")

    # Extract Clipping points
    clipping_points = get_clipping_points(dirpath, point_path)

    # Read and check model
    if not path.exists(model_path):
        RuntimeError("The given directory: %s did not contain the file: model.vtp" % dirpath)

    # Clean surface
    surface = read_polydata(model_path)
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

    #Get a capped and uncapped version of the surface
    open_surface = surface
    capped_surface = capp_surface(surface)

    # Get inlet and outlets
    inlet, outlets = get_centers(open_surface, dirpath)

    # Compute centerline
    centerlines_complete = vmtk_compute_centerlines(inlet, outlets,
                                               centerline_complete_path,
                                               capped_surface, resampling=0.1)
    centerlines_in_order = sort_centerlines(centerlines_complete)

    # Set clipping points in order and make VTK objects
    line = extract_single_line(centerlines_in_order, 0)
    p1, p2, ID1, ID2, vtk_clipping_points, clipping_points = get_vtk_clipping_points(line, clipping_points)
    ID_mid = int((ID1 + ID2) / 2.)

    # Special cases including the opthalmic artery
    eye, clip_ID, centerlines_complete, eyeline = find_ophthalmic_artery(centerlines_complete, clipping_points)

    # Move new line manually
    print("Getting clipped curve.")
    clipped_curve = extract_single_line(centerlines_in_order, 0, startID=ID1, endID=ID2)
    patch_cl = CreateParentArteryPatches(centerlines_in_order, vtk_clipping_points, siphon=True)
    patch_start = extract_single_line(patch_cl, 0)
    patch_ends = []
    n = centerlines_in_order.GetNumberOfCells()
    for i in range(1,n+1):
        patch_ends.append(extract_single_line(patch_cl, i))

    # Find directions to move centerline
    direction = "horizont"
    middle_points, middleIds = get_spline_points(line, beta, direction, vtk_clipping_points)
    dx_p1 = middle_points[0] - p1
    dx_p2 = middle_points[-1] - p2

    print("Moving centerline manually")
    # Move horizontally
    patchline_1 = patch_start
    patchline_2 = patch_ends[0]
    patch_cl_new1 = move_line_horizontally(patchline_1, ID1, ID2, dx_p1, clip=False, side="left")
    patch_cl_new2 = move_line_horizontally(patchline_2, ID1, ID2, dx_p1, clip=False, side="right")
    clipped_part_new = move_line_horizontally(clipped_curve, ID1, ID2, dx_p1, clip=True, eye=eye)
    new_centerline = merge_data([patch_cl_new1, clipped_part_new, patch_cl_new2])
    clipped_part_new = connect_line(clipped_part_new)

    # Move vertically
    direction = "vertical"
    middle_points, middleIds, dx = get_spline_points(new_centerline, alpha, direction, vtk_clipping_points)
    clipped_part_new = move_points_vertically(clipped_part_new, dx)
    new_centerline = merge_data([patch_cl_new1, clipped_part_new, patch_cl_new2])
    new_centerline = connect_line(new_centerline)

    # Extract siphons from old and new centerline
    locator = get_locator(new_centerline)
    ica_siphon = read_polydata(siphon_path)
    endID = ica_siphon.GetNumberOfPoints() - 1
    endID_new  = locator.FindClosestPoint(ica_siphon.GetPoint(endID))
    new_ica_siphon = newsiphon = cut_centerline(new_centerline, endID_new)
    newID1  = locator.FindClosestPoint(p1)
    newID2  = locator.FindClosestPoint(p2)
    newline_sp, newcurv_sp = spline_centerline(new_ica_siphon, get_curv=True, isline=True, nknots=11)
    new_p1 = newline_sp.GetPoints().GetPoint(newID1)
    new_p2 = newline_sp.GetPoints().GetPoint(newID2)
    moved_carotid_siphon = extract_single_line(newline_sp, 0, startID=newID1, endID=newID2)

    # Rename
    siphon = ica_siphon
    moved_siphon = moved_carotid_siphon

    if method in ["maxcurv", "odrline","smooth", "frac"]:
        nknots = 11
        line_sp, curv_sp = spline_centerline(siphon, get_curv=True, isline=True, nknots=nknots)
        curv_sp= resample(curv_sp, line_sp.GetNumberOfPoints())
        newline_sp, newcurv_sp = spline_centerline(newsiphon, get_curv=True, isline=True, nknots=nknots)
        cutcurv = curv_sp[ID1:ID2]
        newcutcurv = newcurv_sp[newID1:newID2]

    if method == "discrete":
        # Smooth line with discrete derivatives
        neigh = 30
        line_d, curv_d =  discrete_geometry(siphon, neigh=neigh)
        newline_d, newcurv_d =  discrete_geometry(newsiphon, neigh=neigh)
        cutcurv_d = curv_d[ID1:ID2]
        newcutcurv_d = newcurv_d[newID1:newID2]


    if method == "MISR":
        # Map MISR values to old and new splined anterior bend
        anterior_bend =extract_single_line( line, 0, startID=ID1, endID=ID2)
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


    if proj:
        print "Computing 2D Angles"
    else:
        print "Computing 3D Angles"

    # Find adjusted clipping points (and tracing points)
    if method == "plane":
        maxP, maxID     = find_furthest_points(dx, siphon)
        newmaxP, newmaxID = find_furthest_points(dx, moved_siphon)

    elif method in ["itplane", "itplane_clip"]:
        maxP, maxID     = find_furthest_points(dx, siphon )
        newmaxP, newmaxID = find_furthest_points(dx, moved_siphon)

        siphon = vmtk_centerline_geometry(siphon, False)

        T1 = get_array("FrenetTangent", siphon, k=3)
        T2 = get_array("FrenetTangent", moved_siphon, k=3)

        p1_1, p1_id = find_closest_point(T1[-1],0, maxID, p2, siphon )
        p2_2, p2_id = find_closest_point(T1[0],maxID, siphon.GetNumberOfPoints(), p1, siphon)

        newp1_1, np1_id = find_closest_point(T2[-1],0, newmaxID , new_p2, moved_siphon)
        newp2_2, np2_id = find_closest_point(T2[0],newmaxID, moved_siphon.GetNumberOfPoints(), new_p1, moved_siphon)

        N1 = get_array("FrenetBinormal", siphon, k=3)[p1_id]
        N2 = get_array("FrenetBinormal", moved_siphon, k=3)[np1_id]

        dP = p1_1 - p2_2
        dnewP = newp1_1 - newp2_2

        normal = np.cross(dP, N1)
        newnormal = np.cross(dnewP, N2)

        maxP, maxID     = find_furthest_points(normal, siphon)
        newmaxP, newmaxID = find_furthest_points(newnormal, moved_siphon)

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

    # Compute angles based on the classic formula for
    # angle between vectors in 3D
    if method == "odrline":
        limits = ["cumulative" , "sd"]
        pA = pB = newpA = newpB = np.zeros(3)
        for limit in limits:
            d1, d2, curvlineold = odr_line(ID1, ID2, line_sp, curv_sp, limit)
            newd1, newd2, curvlinenew = odr_line(newID1, newID2, newline_sp, newcurv_sp, limit)

            deg = find_angle_odr(d1,d2,proj)
            newdeg = find_angle_odr(newd1,newd2,proj)

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


    else:
        if method == "frac":
            n_values = [5] #[3,5,5,7,7]
            l = [2] #[1,2,1,1,2]
            r = [3] ##[2,3,4,6,5]
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

                newdeg, nl1,nl2 = find_angle(newpA, newpB, new_p1,new_p2, proj)
                return newdeg


                old_deg.append(deg)
                new_deg.append(newdeg)

        elif method in ["plane", "itplane", "itplane_clip", "maxcurv", "smooth", "discrete", "maxdist"]:
            frac_vals = [ 2. / 5. ]
            for frac in frac_vals:

                if method == "itplane_clip":
                    IDmid = (p2_id - p1_id)/2.
                    newIDmid = (np2_id - np1_id)/2.
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

                else:
                    IDA = int(maxID*frac)
                    IDB = int(maxID*(1 + (1-frac)))
                    pA = siphon.GetPoints().GetPoint(IDA)
                    pB = siphon.GetPoints().GetPoint(IDB)

                    deg,l1,l2 = find_angle(pA, pB, p1,p2, proj)

                    IDA = int(newmaxID*frac)
                    IDB = int(newmaxID*(1 + (1-frac)))
                    newpA = moved_siphon.GetPoints().GetPoint(IDA)
                    newpB = moved_siphon.GetPoints().GetPoint(IDB)

                    newdeg,nl1,nl2 = find_angle(newpA, newpB, new_p1,new_p2, proj)
    return newdeg


def compute_curvature(dirpath, point_path, name, alpha, beta, method):
    """
    Primary collection of methods for computing curvature of a centerline.
    Five methods are currently implemented:
    1) VMTK - Factor variance (vmtkfactor)
    2) VMTK - Iteration variance (vmtkit)
    3) Discrete derivatives (disc)
    4) Knot free regression splines (knotfree) [Requires MATLAB]
    5) B-splines (spline)

    Args:
        dirpath (str): Path to case folder.
        point_path (str): Filename of clipping point file.
        name (str): Directory where surface models are located.
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        method (str): Method used to compute curvature.

    Returns:
        maxcurv (float): Maximum curvature within the selected siphon.
    """
    # Input filenames
    model_path = path.join(dirpath, name, "model.vtp")
    model_new_surface = path.join(dirpath, name, "new_model_alpha_%s_beta_%s.vtp" % (alpha,beta))

    # Centerliens
    centerline_complete_path = path.join(dirpath, name, "centerline_complete.vtp")
    new_centerlines_path = path.join(dirpath, name, "new_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta))

    # Extract Clipping points
    clipping_points = get_clipping_points(dirpath, point_path)

    # Read and check model
    if not path.exists(model_path):
        RuntimeError("The given directory: %s did not contain the file: model.vtp" % dirpath)

    # Clean surface
    surface = read_polydata(model_path)
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

    #Get a capped and uncapped version of the surface
    open_surface = surface
    capped_surface = capp_surface(surface)

    # Get inlet and outlets
    inlet, outlets = get_centers(open_surface, dirpath)

    # Compute centerline
    centerlines_complete = vmtk_compute_centerlines(inlet, outlets,
                                               centerline_complete_path,
                                               capped_surface, resampling=0.1)

    centerlines_in_order = sort_centerlines(centerlines_complete)

    # Set clipping points in order and make VTK objects
    line = extract_single_line(centerlines_in_order, 0)
    p1, p2, ID1, ID2, vtk_clipping_points, clipping_points = get_vtk_clipping_points(line, clipping_points)
    ID_mid = int((ID1 + ID2) / 2.)

    # Special cases including the opthalmic artery
    eye, clip_ID, centerlines_complete, eyeline = find_ophthalmic_artery(centerlines_complete, clipping_points)

    print("Clipping centerlines.")
    patch_cl = CreateParentArteryPatches(centerlines_in_order, vtk_clipping_points, siphon=True)
    clipped_curve = extract_single_line(centerlines_in_order, 0, startID=ID1, endID=ID2)
    patch_start = extract_single_line(patch_cl, 0)
    patch_ends = []
    n = centerlines_in_order.GetNumberOfCells()
    for i in range(1,n+1):
        patch_ends.append(extract_single_line(patch_cl, i))

    # Find ID of middle pooint:
    direction_tmp = "horizont"
    middle_points, middleIds = get_spline_points(line, beta, direction_tmp,  vtk_clipping_points)
    dx_p1 = middle_points[0] - p1
    dx_p2 = middle_points[-1] - p2
    maxcurv = None

    # Compute new centerline using VMTK
    if "vmtk" in method:
        new_centerline_vmtk = make_centerline(model_new_surface, new_centerlines_path, smooth=False, resampling=False)
        centerlines_in_order = sort_centerlines(new_centerline_vmtk)
        new_line_vmtk = extract_single_line(centerlines_in_order, 0)

        # Find new boundaries
        locator = get_locator(new_line_vmtk)
        ID1_new = locator.FindClosestPoint(p1)
        ID2_new = locator.FindClosestPoint(p2)

        # 1) VMTK - Factor variance
        if method == "vmtkfactor":
            factor = 0.5
            line_fac = vmtk_centerline_geometry(new_line_vmtk, smooth=True, iterations=100, factor=factor)
            curv_fac = get_array("Curvature", line_fac)
            curv_fac= gauss(curv_fac, 5)
            maxcurv = max(curv_fac[ID1_new+10:ID2_new-10])[0]


        # 2) VMTK - Iteration variance
        elif method == "vmtkit":
            it = 150
            line_it = vmtk_centerline_geometry(new_line_vmtk, smooth=True, iterations=it, factor=1.0)
            curv_it = get_array("Curvature", line_it)
            curv_it= gauss(curv_it, 5)
            maxcurv = max(curv_it[ID1_new+10:ID2_new-10])[0]

    else:
        # Compute new centerline by manual displacement
        print("Moving centerline manually")
        p_1 = patch_start
        p_2 = patch_ends[0]
        patch_cl_new1 = move_line_horizontally(p_1, ID1, ID2, dx_p1, clip=False, side="left")
        patch_cl_new2 = move_line_horizontally(p_2, ID1, ID2, dx_p1, clip=False, side="right")
        clipped_part_new = move_line_horizontally(clipped_curve, ID1, ID2, dx_p1, clip=True, eye=eye)
        new_centerline = merge_data([patch_cl_new1, clipped_part_new, patch_cl_new2])
        clipped_part_new = connect_line(clipped_part_new)
        direction = "vertical"
        middle_points, middleIds, dx = get_spline_points(new_centerline, alpha, direction, vtk_clipping_points)
        clipped_part_new = move_points_vertically(clipped_part_new, dx)
        new_centerline = merge_data([patch_cl_new1, clipped_part_new, patch_cl_new2])
        new_centerline = connect_line(new_centerline)

        # Find new boundaries of siphon
        locator = get_locator(new_centerline)
        ID1_new = locator.FindClosestPoint(p1)
        ID2_new = locator.FindClosestPoint(p2)

        # 3) Discrete derivatives
        if method == "disc":
            neigh = 20
            line_di, curv_di = discrete_geometry(new_centerline, neigh=neigh)
            filtercurv = gauss(curv_di,5)
            maxcurv = max(filtercurv[ID1_new+10:ID2_new-10])

        # 5) Knot free regression splines
        if method == "knotfree":
            clfile = write_centerline(new_centerline)
            inits = 20
            curv_kf = spline_matlab_noload(".", clfile, init_knots=inits, order=float(5))
            curv_kf = gauss(curv_kf,5)
            maxcurv = max(curv_kf[ID1_new+10:ID2_new-10])

        # 6) Splines
        if method == "spline":
            nknots = 50
            line_sp, curv_sp = spline_centerline(new_centerline, get_curv=True, isline=True, nknots=nknots)
            curv_sp = gauss(curv_sp,5)
            maxcurv = max(curv_sp[ID1_new+10:ID2_new-10])

    return maxcurv


def odr_line(ID1, ID2, line, curvature, limit):
    """
    Computes the othogonal distance regression
    of points along the centerline selected from
    1) All points until a cumulative limit is reached
    or
    2) The first 11 points and all points fulfilling curvature
    less than the mean plus 1.96 x SD

    Args:
        ID1 (int): ID of first clipping point.
        ID2 (int): ID of second clipping point.
        line (vtkPolyData): Centerline data.
        curvature (ndarray): Array of curvature values.
        limit (ndarray): Method used as limit

    Returns:
        d1 (ndarray): Direction vector from first clipping point.
        d2 (ndarray): Direction vector from second clipping point.
        curvlines (vtkPolyData): Centerline object with corresponding curvature values.
    """
    lim = len(curvature)-1

    if limit == "cumulative":
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
    return d1, d2, curvlines

def cut_centerline(line, endID):
    """
    Clips the centerline at the end point
    given by some ID.

    Args:
        line (vtkPolyData): Centerline data.
        endID (int): ID of point to cut.

    Returns:
        line (vtkPolyData): Clipped line.
    """
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


def find_angle(pA, pB, p1,p2, proj):
    """
    Compute the angle between two vectors
    a = pA - p1 and b = pB - p2
    using the classical formula.

    Args:
        pA (ndarray): Point along the centerline.
        pB (ndarray): Point along the centerline.
        p1 (ndarray): Point along the centerline.
        p2 (ndarray): Point along the centerline.
        proj (bool): True / False for 2D / 3D angle.

    Returns:
        deg (float): Angle between vectors.
        P1A (ndarray): First vector.
        P2B (ndarraty): Second vector.
    """
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
    """
    Compute the angle between two vectors
    d1 and d2 using the classical formula.
    Used for the ODR-method, spesifically.

    Args:
        d1 (ndarray): First vector
        d2 (ndarray): Second vector
        proj (bool): True / False for 2D / 3D angle.

    Returns:
        deg (float): Angle between vectors.
    """
    if d1.dot(d2) > 0:
        d1 = -d1
    if proj:
        d1[0] = 0
        d2[0] = 0

    costheta = (d1.dot(-d2)) / (la.norm(d1)*la.norm(-d2))
    angle = np.arccos(costheta)
    deg = (angle * 180 / np.pi)
    return deg, d1, d2

def save_angle_or_curvature(values, case, param):
    """
    Save values of curvature / angle stored in a
    n x n matrix.

    Args:
        values (ndarray): n x n matrix containing values.
        case (str): Name of case.
        param (str): Name of parameter stored.
    """
    mat = np.matrix(values)
    with open('new_%s_%s.txt' % (param,case),'wb') as f:
        for line in mat:
            np.savetxt(f, line, fmt='%.3f')

def spline_matlab_noload(path, filename, init_knots, order):
    """
    Perform Knot-free regresion spline on input centerline
    extracted from input filename.
    Excludes loading of matlab engine.

    Args:
        path (str): Path to centerline-text location.
        filename (str): Filename of text where centerlines is stored.
        init_knots (int): Number of initial knots.
        order (int): Order of spline.

    Returns:
        curv_p (ndarray): Array of curvature values.
    """

    print("Computing knot free regression spline")
    order = float(order)
    init_knots = float(init_knots)
    curv_m = mlab.CenterlineCharacterization(path, filename, init_knots, order, nargout=1)

    n = len(curv_m)
    curv_p = np.zeros(n)

    for i in range(n):
        curv_p[i] = curv_m[i][0]

    return curv_p


def main(basedir, case, kappa, theta, alpha, beta, method_curv, method_angle, n=50):
    """
    Initilization for computing curvature and angle.
    Values are either printed to terminal or stored in a (n x n) matrix.

    Args:
        basedir (str): Location of case folders.
        case (str): Name of case.
        kappa (bool): True to compute curvature.
        theta (bool): True to compute angle.
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        method_curv (str): Method used to compute curvature.
        method_angle (str): Method used to compute angle.
        n (int): Determines matrix size when computing multiple values.
    """
    name = "surface"
    point_path = "carotid_siphon_points.particles"

    # One or more cases
    if case is not None:
        folders = [case]
    else:
        folders = sorted([folder for folder in listdir(basedir) if folder[:2] in ["P0"]])

    # Movement in one or multiple directions
    if alpha is not None:
        alphas = [alpha]
        betas = [beta]
    else:
        ab_bound = np.loadtxt("alphabeta_bound.txt")
        max_curv_values = np.zeros((n,n))
        angle_values = np.zeros((n,n))
        k = 0

    # Iterate through cases and compute quantities
    for folder in folders:
        print("Working on case " + folder)
        casedir = path.join(basedir,folder)

        if alpha is None:
            amin,amax,bmin,bmax = ab_bound[k][0], ab_bound[k][1], ab_bound[k][2], ab_bound[k][3]
            alphas = np.linspace(amin,amax, n)
            betas = np.linspace(bmin,bmax,n)

        for i,alpha in enumerate(alphas):
            for j,beta in enumerate(betas):
                # Compute curvature (kappa) or / and angle (theta)
                if kappa:
                    maxcurv = compute_curvature(casedir, point_path, name,  alpha,beta, method_curv )
                if theta:
                    angle = compute_angle(casedir, point_path, name,alpha, beta, method_angle)

                if len(alphas) > 1:
                    if kappa:
                        max_curv_values[i,j] = maxcurv
                    if theta:
                        angle_values[i,j] = angle
                else:
                    if kappa:
                        print("Curvature = %.3f" % maxcurv)
                    if theta:
                        print("Angle = %.3f" % angle)

        if len(alphas) > 1:
            if kappa:
                save_angle_or_curvature(max_curv_values, folder, "curvature")
            if theta:
                save_angle_or_curvature(angle_values, folder, "angle")
            k+=1




if  __name__ == "__main__":
    basedir, case, kappa, theta, alpha, beta, method_curv, method_angle = read_command_line()
    if method_curv == "knotfree":
        mlab = matlab.engine.start_matlab()
    main(basedir, case, kappa, theta, alpha, beta, method_curv, method_angle)

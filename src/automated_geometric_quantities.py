from argparse import ArgumentParser, RawDescriptionHelpFormatter
from scipy.signal import argrelextrema
from scipy.ndimage.filters import gaussian_filter as gauss
import operator
# Local import
from common import *


def compute_angle(input_filepath, alpha, beta, method, new_centerlines, proj=False):
    """
    Primary collection of methods for computing the angle of a vessel bend.
    Three main methods are currently implemented:
    1) ODR methods: odrline
    2) Tracing point methods: maxcurv, smooth, discrete, frac, MISR
    3) Relative tracing point methods: plane, itplane, itplane_clip

    Args:
        input_filepath (str): Path to case folder.
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        method (str): Method used to compute angle.
        proj (bool): True / False for computing 2D / 3D angle.
        new_centerlines (vtkPolyData): New centerline.

    Returns:
        newdeg (float): New angle of a vessel bend from a manipulated centerline.
    """
    # Case name

    print("Computing angle")
    base_path = get_path_names(input_filepath)

    # Find endID from landmarking
    centerline_path = base_path + "_centerline.vtp"

    # Extract Clipping points
    point_path = base_path + "_anterior_bend.particles"
    if not path.exists(point_path):
        RuntimeError("The given .particles file: %s does not exist!" % point_path)
    region_points = np.loadtxt(point_path)
    p1 = region_points[0]
    p2 = region_points[1]

    # Extract old centerline
    if not path.exists(centerline_path):
        RuntimeError("The given .vtp file: %s does not exist!" % centerline_path)
    centerlines = read_polydata(centerline_path)

    if new_centerlines is None:
        print("Maniuplating centerline manually")
        centerlines, new_centerlines = get_new_centerlines(centerlines, region_points, alpha, beta, p1, p2)

    # Get new siphon and prepare
    id1, id2, moved_id1, moved_id2, moved_p1, moved_p2 = get_moved_siphon(new_centerlines, centerlines, p1, p2)

    # Get anterior bend only
    siphon = extract_single_line(centerlines, 1, startID=id1, endID=id2)
    moved_siphon = extract_single_line(new_centerlines, 0, startID=moved_id1, endID=moved_id2)
    id1, id2 = 0, siphon.GetNumberOfPoints() - 1
    moved_id1, moved_id2 = 0, moved_siphon.GetNumberOfPoints() - 1

    if method in ["maxcurv", "odrline", "smooth", "frac"]:
        nknots = 11
        siphon_splined, siphon_curv = spline_centerline(siphon, get_curv=True, isline=True, nknots=nknots,
                                                        get_misr=False)

        moved_siphon_splined, moved_siphon_curv = spline_centerline(moved_siphon, get_curv=True, isline=True,
                                                                    nknots=nknots, get_misr=False)
        siphon_curv = resample(siphon_curv, siphon_splined.GetNumberOfPoints())
        cutcurv = siphon_curv[id1:id2]
        newcutcurv = moved_siphon_curv[moved_id1:moved_id2]

    if method == "discrete":
        # Smooth line with discrete derivatives
        neigh = 30
        line_d, curv_d = discrete_geometry(siphon, neigh=neigh)
        newline_d, newcurv_d = discrete_geometry(moved_siphon, neigh=neigh)
        cutcurv_d = curv_d[id1:id2]
        newcutcurv_d = newcurv_d[moved_id1:moved_id2]

    if method == "MISR":
        # Map MISR values to old and new splined anterior bend
        anterior_bend = extract_single_line(centerlines, 0, startID=id1, endID=id2)
        m = anterior_bend.GetNumberOfPoints()
        m1 = moved_siphon.GetNumberOfPoints()
        misr_array = get_vtk_array(radiusArrayName, 1, m)
        newmisr_array = get_vtk_array(radiusArrayName, 1, m1)
        misr_list = []
        for i in range(m):
            misr = anterior_bend.GetPointData().GetArray(radiusArrayName).GetTuple(i)
            misr_list.append(misr[0])
            misr_array.SetTuple(i, misr)

        misr_list = resample(misr_list, m1)
        for i in range(m1):
            newmisr_array.SetTuple(i, (misr_list[i],))

        siphon.GetPointData().AddArray(misr_array)
        moved_siphon.GetPointData().AddArray(newmisr_array)

    if proj:
        print("Computing 2D Angles")
    else:
        print("Computing 3D Angles")

    # Get direction to the point furthest away (dx)
    direction = "vertical"
    clipping_points_vtk = vtk.vtkPoints()
    for point in np.asarray(region_points):
        clipping_points_vtk.InsertNextPoint(point)
    middle_points, middle_ids, dx = get_spline_points(extract_single_line(centerlines, 0), 0.1, direction,
                                                      clipping_points_vtk)

    # Find adjusted clipping points (and tracing points)
    if method == "plane":
        max_p, max_id = find_furthest_points(dx, siphon)
        newmax_p, newmax_id = find_furthest_points(dx, moved_siphon)

    elif method in ["itplane", "itplane_clip"]:
        max_p, max_id = find_furthest_points(dx, siphon)
        newmax_p, newmax_id = find_furthest_points(dx, moved_siphon)

        siphon = vmtk_centerline_geometry(siphon, False)

        T1 = get_array("FrenetTangent", siphon, k=3)
        T2 = get_array("FrenetTangent", moved_siphon, k=3)

        p1_1, p1_id = find_closest_point(T1[-1], 0, max_id, p2, siphon)
        p2_2, p2_id = find_closest_point(T1[0], max_id, siphon.GetNumberOfPoints(), p1, siphon)

        newp1, np1_id = find_closest_point(T2[-1], 0, newmax_id, moved_p2, moved_siphon)
        newp2, np2_id = find_closest_point(T2[0], newmax_id,
                                           moved_siphon.GetNumberOfPoints(), moved_p1,
                                           moved_siphon)

        n1 = get_array("FrenetBinormal", siphon, k=3)[p1_id]
        n2 = get_array("FrenetBinormal", moved_siphon, k=3)[np1_id]

        dp = p1_1 - p2_2
        dnewp = newp1 - newp2

        normal = np.cross(dp, n1)
        newnormal = np.cross(dnewp, n2)

        max_p, max_id = find_furthest_points(normal, siphon)
        newmax_p, newmax_id = find_furthest_points(newnormal, moved_siphon)

    elif method == "maxcurv":
        max_id, v = max(enumerate(cutcurv), key=operator.itemgetter(1))
        newmax_id, v = max(enumerate(newcutcurv), key=operator.itemgetter(1))
        max_p = siphon_splined.GetPoint(id1 + max_id)
        newmax_p = moved_siphon.GetPoint(moved_id1 + newmax_id)

    elif method == "smooth":
        allmaxcurv = argrelextrema(cutcurv, np.greater)[0]
        allnewmaxcurv = argrelextrema(newcutcurv, np.greater)[0]

        tmpcurv = cutcurv
        while len(allmaxcurv) > 2:
            tmpcurv = gauss(tmpcurv, 2)
            allmaxcurv = argrelextrema(tmpcurv, np.greater)[0]

        tmpnewcurv = newcutcurv
        while len(allnewmaxcurv) > 2:
            tmpnewcurv = gauss(tmpnewcurv, 2)
            allnewmaxcurv = argrelextrema(tmpnewcurv, np.greater)[0]

        max_id = allmaxcurv[0]
        newmax_id = allnewmaxcurv[0]

    elif method == "discrete":
        max_id, v = max(enumerate(cutcurv_d), key=operator.itemgetter(1))
        newmax_id, v = max(enumerate(newcutcurv_d), key=operator.itemgetter(1))

    elif method == "maxdist":
        norm_p1 = [la.norm(np.array(p1) - np.array(siphon.GetPoint(i))) for i in range(siphon.GetNumberOfPoints())]
        norm_p2 = [la.norm(np.array(p2) - np.array(siphon.GetPoint(i))) for i in
                   range(siphon.GetNumberOfPoints() - 1, -1, -1)]
        max_id = 0
        max_dist = 0
        for i, n1 in enumerate(norm_p1):
            for j, n2 in enumerate(norm_p2):
                dist = n1 ** 2 + n2 ** 2
                if dist > max_dist:
                    max_dist = dist
                    max_id = i

        newnorm_p1 = [la.norm(np.array(moved_p1) - np.array(moved_siphon.GetPoint(i))) for i in
                      range(moved_siphon.GetNumberOfPoints())]
        newnorm_p2 = [la.norm(np.array(moved_p2) - np.array(moved_siphon.GetPoint(i))) for i in
                      range(moved_siphon.GetNumberOfPoints() - 1, -1, -1)]
        newmax_id = 0
        new_max_dist = 0
        for i, n1 in enumerate(newnorm_p1):
            for j, n2 in enumerate(newnorm_p2):
                dist = n1 ** 2 + n2 ** 2
                if dist > new_max_dist:
                    new_max_dist = dist
                    newmax_id = i

    # Compute angles based on the classic formula for
    # angle between vectors in 3D
    if method == "odrline":
        limits = ["cumulative", "sd"]
        for limit in limits:
            d1, d2, curvlineold = odr_line(id1, id2, siphon_splined, siphon_curv, limit)
            newd1, newd2, curvlinenew = odr_line(moved_id1, moved_id2, moved_siphon, moved_siphon_curv, limit)

            deg = find_angle_odr(d1, d2, proj)
            newdeg = find_angle_odr(newd1, newd2, proj)

    elif method == "MISR":
        multiplier = 1.5
        N1 = siphon.GetNumberOfPoints()
        N2 = moved_siphon.GetNumberOfPoints()
        rad1 = siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(0)
        rad2 = siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(N1 - 1)
        newrad1 = moved_siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(0)
        newrad2 = moved_siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(N2 - 1)

        pA, rA = move_past_sphere(siphon, p1, rad1, 0, step=1, stop=N1 - 1, X=multiplier)
        pB, rB = move_past_sphere(siphon, p2, rad2, N1 - 1, step=-1, stop=0, X=multiplier)
        newpA, rA = move_past_sphere(moved_siphon, moved_p1, newrad1, 0, step=1, stop=N2 - 1, X=multiplier)
        newpB, rB = move_past_sphere(moved_siphon, moved_p2, newrad2, N2 - 1, step=-1, stop=0, X=multiplier)

        deg, l1, l2 = find_angle(pA, pB, p1, p2, proj)
        newdeg, nl1, nl2 = find_angle(newpA, newpB, moved_p1, moved_p2, proj)

    else:
        if method == "frac":
            n_values = [5]
            l = [2]
            r = [3]
            i = 0
            dX = 1. / n_values[i]
            IDA = int(id1 + (id2 - id1) * l[i] * dX)
            IDB = int(id1 + (id2 - id1) * r[i] * dX)
            pA = siphon_splined.GetPoints().GetPoint(IDA)
            pB = siphon_splined.GetPoints().GetPoint(IDB)

            IDA = int(moved_id1 + (moved_id2 - moved_id1) * l[i] * dX)
            IDB = int(moved_id1 + (moved_id2 - moved_id1) * r[i] * dX)
            newpA = moved_siphon.GetPoints().GetPoint(IDA)
            newpB = moved_siphon.GetPoints().GetPoint(IDB)

            deg, l1, l2 = find_angle(pA, pB, p1, p2, proj)
            newdeg, nl1, nl2 = find_angle(newpA, newpB, moved_p1, moved_p2, proj)

        elif method in ["plane", "itplane", "itplane_clip", "maxcurv", "smooth",
                        "discrete", "maxdist"]:
            frac = 4. / 5.
            if method == "itplane_clip":
                IDmid = (p2_id - p1_id) / 2.
                newIDmid = (np2_id - np1_id) / 2.
                if max_id > IDmid:
                    IDA = int((max_id - p1_id) * frac)
                    IDB = int((max_id - p1_id) * (1 + (1 - frac)))
                    pA = siphon.GetPoints().GetPoint(IDA + p1_id)
                    pB = siphon.GetPoints().GetPoint(IDB + p1_id)
                else:
                    IDB = int((p2_id - max_id) * (1 + (1 - frac)))
                    IDA = int((p2_id - max_id) * frac)
                    pA = siphon.GetPoints().GetPoint(IDA)
                    pB = siphon.GetPoints().GetPoint(IDB)

                if newmax_id > newIDmid:
                    IDA = int((newmax_id - np1_id) * frac)
                    IDB = int((newmax_id - np1_id) * (1 + (1 - frac)))
                    newpA = moved_siphon.GetPoints().GetPoint(IDA + np1_id)
                    newpB = moved_siphon.GetPoints().GetPoint(IDB + np1_id)
                else:
                    IDA = int((np2_id - newmax_id) * frac)
                    IDB = int((np2_id - newmax_id) * (1 + (1 - frac)))
                    newpA = moved_siphon.GetPoints().GetPoint(IDA)
                    newpB = moved_siphon.GetPoints().GetPoint(IDB)

                deg, l1, l2 = find_angle(pA, pB, p1, p2, proj)
                newdeg, nl1, nl2 = find_angle(newpA, newpB, moved_p1, moved_p2, proj)

            else:
                IDA = int(max_id * frac)
                IDB = int(max_id * (1 + (1 - frac)))
                pA = siphon.GetPoints().GetPoint(IDA)
                pB = siphon.GetPoints().GetPoint(IDB)

                IDA = int(newmax_id * frac)
                IDB = int(newmax_id * (1 + (1 - frac)))
                newpA = moved_siphon.GetPoints().GetPoint(IDA)
                newpB = moved_siphon.GetPoints().GetPoint(IDB)

                deg, l1, l2 = find_angle(pA, pB, p1, p2, proj)
                newdeg, nl1, nl2 = find_angle(newpA, newpB, moved_p1, moved_p2, proj)

    return newdeg, deg


def compute_curvature(input_filepath, alpha, beta, method, new_centerlines):
    """
    Primary collection of methods for computing curvature of a centerline.
    Five methods are currently implemented:
    1) VMTK - Factor variance (vmtkfactor)
    2) VMTK - Iteration variance (vmtkit)
    3) Discrete derivatives (disc)
    4) B-splines (spline)

    Args:
        input_filepath (str): Path to case folder.
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        method (str): Method used to compute angle.
        new_centerlines (vtkPolyData): New centerline.

    Returns:
        maxcurv (float): Maximum curvature within the selected siphon.
    """
    # Input filename
    base_path = get_path_names(input_filepath)

    # Centerline filename
    centerline_path = base_path + "_centerline.vtp"

    # Extract Clipping points
    point_path = base_path + "_anterior_bend.particles"
    if not path.exists(point_path):
        RuntimeError("The given .particles file: %s does not exist!" % point_path)
    region_points = np.loadtxt(point_path)
    p1 = region_points[0]
    p2 = region_points[1]

    # Extract old centerline
    if not path.exists(centerline_path):
        RuntimeError("The given .vtp file: %s does not exist!" % centerline_path)
    centerlines = read_polydata(centerline_path)

    if new_centerlines is None:
        print("Maniuplating centerline manually")
        centerlines, new_centerlines = get_new_centerlines(centerlines, region_points, alpha, beta, p1, p2)

    # Compute new centerline using VMTK
    new_centerline = extract_single_line(new_centerlines, 0)
    centerline = extract_single_line(centerlines, 0)
    new_locator = get_locator(new_centerline)
    old_locator = get_locator(centerline)
    id1 = old_locator.FindClosestPoint(p1)
    id2 = old_locator.FindClosestPoint(p2)
    id1_new = new_locator.FindClosestPoint(p1)
    id2_new = new_locator.FindClosestPoint(p2)

    # 1) VMTK - Factor variance
    if method == "vmtkfactor":
        factor = 0.5
        line_fac = vmtk_centerline_geometry(new_centerline, smooth=True, iterations=100, factor=factor)
        curv_fac = get_array("Curvature", line_fac)
        new_curvature = gauss(curv_fac, 5)

        line_fac = vmtk_centerline_geometry(centerlines, smooth=True, iterations=100, factor=factor)
        curv_fac = get_array("Curvature", line_fac)
        curvature = gauss(curv_fac, 5)

    # 2) VMTK - Iteration variance
    elif method == "vmtkit":
        it = 150
        line_it = vmtk_centerline_geometry(new_centerline, smooth=True, iterations=it, factor=1.0)
        curv_it = get_array("Curvature", line_it)
        new_curvature = gauss(curv_it, 5)

        line_it = vmtk_centerline_geometry(centerlines, smooth=True, iterations=it, factor=1.0)
        curv_it = get_array("Curvature", line_it)
        curvature = gauss(curv_it, 5)

    # 3) Splines
    elif method == "spline":
        nknots = 50
        siphon_splined, siphon_curv = spline_centerline(new_centerline, get_curv=True,
                                                        isline=True, nknots=nknots)
        new_curvature = gauss(siphon_curv, 5)

        siphon_splined, siphon_curv = spline_centerline(centerlines, get_curv=True,
                                                        isline=True, nknots=nknots)
        curvature = gauss(siphon_curv, 5)

    # 4) Default: Discrete derivatives
    elif method == "disc":
        neigh = 20
        line_di, curv_di = discrete_geometry(new_centerline, neigh=neigh)
        new_curvature = gauss(curv_di, 5)

        line_di, curv_di = discrete_geometry(centerlines, neigh=neigh)
        curvature = gauss(curv_di, 5)

    old_maxcurv = max(curvature[id1 + 10:id2 - 10])
    new_maxcurv = max(new_curvature[id1_new + 10:id2_new - 10])

    return old_maxcurv, new_maxcurv


def get_new_centerlines(centerlines, region_points, alpha, beta, p1, p2):
    centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        find_region_of_interest_and_diverging_centerlines(centerlines, region_points)
    new_centerlines = centerlines
    diverging_id = None if len(diverging_ids) == 0 else diverging_ids[0]
    if beta != 0.0:
        direction = "horizont"
        middle_points, middle_ids = get_spline_points(new_centerlines, beta, direction, region_points_vtk)
        dx_p1 = middle_points[0] - p1
        new_centerlines = move_centerlines(new_centerlines, dx_p1, p1, p2, diverging_id,
                                           diverging_centerlines, direction)
    if alpha != 0.0:
        direction = "vertical"
        middle_points, middle_ids, dx = get_spline_points(new_centerlines, alpha, direction,
                                                          region_points_vtk)
        new_centerlines = move_centerlines(new_centerlines, dx, p1, p2, diverging_id, diverging_centerlines, direction)
    return centerlines, new_centerlines


def odr_line(id1, id2, line, curvature, limit):
    """
    Computes the othogonal distance regression
    of points along the centerline selected from
    1) All points until a cumulative limit is reached
    or
    2) The first 11 points and all points fulfilling curvature
    less than the mean plus 1.96 x SD

    Args:
        id1 (int): ID of first clipping point.
        id2 (int): ID of second clipping point.
        line (vtkPolyData): Centerline data.
        curvature (ndarray): Array of curvature values.
        limit (ndarray): Method used as limit

    Returns:
        d1 (ndarray): Direction vector from first clipping point.
        d2 (ndarray): Direction vector from second clipping point.
        curvlines (vtkPolyData): Centerline object with corresponding curvature values.
    """
    lim = len(curvature) - 1

    if limit == "cumulative":
        max_cum = 10
        id1_up = id1 + 1
        id1_down = id1 - 1
        id2_up = id2 - 1
        id2_down = id2 + 1
        while sum(curvature[id1:id1_up + 1]) < max_cum and id1_up < lim:
            id1_up += 1
        while sum(curvature[id1_down:id1 + 1]) < max_cum and id1_down > 0:
            id1_down -= 1
        while sum(curvature[id2_up:id2 + 1]) < max_cum and id2_up > 0:
            id2_up -= 1
        while sum(curvature[id2:id2_down + 1]) < max_cum and id2_down < lim:
            id2_down += 1
    else:
        id1_up = id1 + 5
        id1_down = id1 - 5
        id2_up = id2 - 5
        id2_down = id2 + 5

        mean1 = sum(curvature[id1_down:id1_up + 1]) / 11.
        mean2 = sum(curvature[id2_up:id2_down + 1]) / 11.

        sd_1 = np.sqrt(sum((curvature[id1_down:id1_up + 1] - mean1) ** 2) / 10)
        sd_2 = np.sqrt(sum((curvature[id2_up:id1_down + 1] - mean2) ** 2) / 10)
        tol1 = mean1 + sd_1 * 1.96
        tol2 = mean2 + sd_2 * 1.96

        while curvature[id1_up] < tol1 and id1_up < lim:
            id1_up += 1
        while curvature[id1_down] < tol1 and id1_down > 0:
            id1_down -= 1
        while curvature[id2_up] < tol2 and id2_up > 0:
            id2_up -= 1
        while curvature[id2_down] < tol2 and id2_down < lim:
            id2_down += 1

    p1s = []
    for i in range(id1_down, id1_up + 1):
        p1s.append(line.GetPoint(i))

    p2s = []
    for i in range(id2_up, id2_down + 1):
        p2s.append(line.GetPoint(i))

    # Arrange points in matrix
    x_1 = np.array([list(p) for p in p1s])
    x_2 = np.array([list(p) for p in p2s])

    # Find mean of points
    avg1 = np.array([np.mean(x_1[:, 0]), np.mean(x_1[:, 1]), np.mean(x_1[:, 2])])
    avg2 = np.array([np.mean(x_2[:, 0]), np.mean(x_2[:, 1]), np.mean(x_2[:, 2])])

    # Subtract the mean from all points
    dx_1 = x_1 - np.array([avg1] * len(x_1))
    dx_2 = x_2 - np.array([avg2] * len(x_2))

    # Find SVD
    _, _, v_1 = la.svd(dx_1)
    _, _, v_2 = la.svd(dx_2)
    # Find direction vector
    d1 = v_1[0]
    d2 = v_2[0]

    # Parametric equation P = p0 + t*d
    # Make lines with curv
    # Create edges between new_centerline points
    curv_lines_split = []
    points = [p1s, p2s]
    for k, p in enumerate(points):
        pts = vtk.vtkPoints()
        for i in range(len(p)):
            pts.InsertNextPoint(p[i])

        lines = vtk.vtkCellArray()
        for i in range(len(p) - 2):
            newline = vtk.vtkLine()
            newline.GetPointIds().SetId(0, i)
            newline.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(newline)

        line_ = vtk.vtkPolyData()
        line_.SetPoints(pts)
        line_.SetLines(lines)

        m = line_.GetNumberOfPoints()
        curv_array = get_vtk_array("Curvature", 1, m)
        if k == 0:
            for i in range(id1_up + 1 - id1_down):
                curv_array.SetTuple(i, [curvature[id1_down + i]])
        else:
            for i in range(id2_down + 1 - id2_up):
                curv_array.SetTuple(i, [curvature[id2_up + i]])

        line_.GetPointData().AddArray(curv_array)

        curv_lines_split.append(line_)

    curvlines = merge_data(curv_lines_split)

    return d1, d2, curvlines


def get_moved_siphon(new_centerlines, centerlines, p1, p2):
    """
    Extracts new siphon from new centerline
    and clipping point information.

    Args:
        new_centerlines (vtkPolyData): Centerline data.
        centerlines (vtkPolyData): Initial centerlines.
        p1 (ndarray): First clipping point.
        p2 (ndarray): Second clipping point.

    Returns:
        moved_id1 (int): New ID of first clipping point.
    Returns:
        moved_id1 (int): New ID of secon clipping point.
    Returns:
        moved_p1 (ndarray): New position of first clipping point.
    Returns:
        moved_p2 (ndarray): New position of second ipping point.
    Returns:
        moved_siphon (vtkPolyData): Splined siphon centerline.
    Returns:
        moved_siphon_curv (ndarray): Curvature array along siphon.
    """
    # Extract new siphon and prepare
    new_locator = get_locator(extract_single_line(new_centerlines, 0))
    old_locator = get_locator(extract_single_line(centerlines, 0))
    id1 = old_locator.FindClosestPoint(p1)
    id2 = old_locator.FindClosestPoint(p2)

    moved_id1 = new_locator.FindClosestPoint(p1)
    moved_id2 = new_locator.FindClosestPoint(p2)
    moved_p1 = new_centerlines.GetPoints().GetPoint(moved_id1)
    moved_p2 = new_centerlines.GetPoints().GetPoint(moved_id2)
    return id1, id2, moved_id1, moved_id2, moved_p1, moved_p2


def find_angle(pA, pB, p1, p2, proj):
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
        vector_a (ndarray): First vector.
        vector_b (ndarraty): Second vector.
    """
    if not proj:
        vector_a = np.array([pA[0] - p1[0], pA[1] - p1[1], pA[2] - p1[2]])
        vector_b = np.array([pB[0] - p2[0], pB[1] - p2[1], pB[2] - p2[2]])
    else:
        vector_a = np.array([0, pA[1] - p1[1], pA[2] - p1[2]])
        vector_b = np.array([0, pB[1] - p2[1], pB[2] - p2[2]])
    costheta = (vector_a.dot(vector_b)) / (la.norm(vector_a) * la.norm(vector_b))
    angle = np.arccos(costheta)
    deg = (angle * 180 / np.pi)

    return deg, vector_a, vector_b


def find_angle_odr(d1, d2, proj):
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

    costheta = (d1.dot(-d2)) / (la.norm(d1) * la.norm(-d2))
    angle = np.arccos(costheta)
    deg = (angle * 180 / np.pi)

    return deg, d1, d2


def save_angle_or_curvature(values, param):
    """
    Save values of curvature / angle stored in a
    n x n matrix.

    Args:
        values (ndarray): n x n matrix containing values.
        param (str): Name of parameter stored.
    """
    mat = np.matrix(values)
    with open('new_%s_%s.txt' % param, 'wb') as f:
        for line in mat:
            np.savetxt(f, line, fmt='%.3f')


def initialize(input_filepath, kappa, theta, alpha, beta, boundary,
               method_curv, method_angle, grid_size=50, proj=False):
    """
    Initilization for computing curvature and angle.
    Values are either printed to terminal or stored in a (n x n) matrix.

    Args:
        input_filepath (str): Location of surface model.
        kappa (bool): True to compute curvature.
        theta (bool): True to compute angle.
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        boundary (list): Bounds of grid for computing quantities.
        method_curv (str): Method used to compute curvature.
        method_angle (str): Method used to compute angle.
        grid_size (int): Determines matrix size when computing multiple values.
        proj (bool): Projects angle into 2d plane if True.
    """
    # New centerlines
    new_centerlines_path = input_filepath[:-4] + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines = read_polydata(new_centerlines_path)

    # Movement in one or multiple directions
    if alpha is not None:
        alphas = [alpha]
        betas = [beta]
    else:
        n = grid_size
        max_curv_values = np.zeros((n, n))
        angle_values = np.zeros((n, n))
        amin, amax, bmin, bmax = boundary[0], boundary[1], boundary[2], boundary[3]
        alphas = np.linspace(amin, amax, n)
        betas = np.linspace(bmin, bmax, n)

    for i, alpha in enumerate(alphas):
        for j, beta in enumerate(betas):
            # Compute curvature (kappa) or / and angle (theta)
            if kappa:
                maxcurv = compute_curvature(input_filepath, alpha, beta, method_curv, new_centerlines)

            if theta:
                angle, oldangle = compute_angle(input_filepath, alpha, beta, method_angle, new_centerlines, proj)

            if len(alphas) > 1:
                if kappa:
                    max_curv_values[i, j] = maxcurv
                if theta:
                    angle_values[i, j] = angle
            else:
                if kappa:
                    print("Curvature = %.3f" % maxcurv)
                if theta:
                    print("Angle = %.3f" % angle)

    if len(alphas) > 1:
        if kappa:
            save_angle_or_curvature(max_curv_values, "curvature")
        if theta:
            save_angle_or_curvature(angle_values, "angle")


def read_command_line():
    """
    Read arguments from commandline
    """
    description = "Computes selected quantities: anterior bend angle or maximum " + \
                  "anterior bend curvature for a given geometry, using the centerline " + \
                  "as a proxy. Possibility for computing several discrete values if "

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required named arguments')

    # Required arguments
    required.add_argument('-i', '--ifile', type=str, required=True,
                          help="Path to the folder with all the cases")
    # Optional arguments
    parser.add_argument('-k', '--curvature', type=bool, default=False,
                        help="Compute curvature variation", metavar="curvature")
    parser.add_argument('-t', '--angle', type=bool, default=False,
                        help="Compute angle variation", metavar="angle")
    parser.add_argument("-a", "--alpha", type=float, default=None,
                        help="Compression factor in vertical direction, ranging from -1.0 to 1.0")
    parser.add_argument("-b", "--beta", type=float, default=None,
                        help="Compression factor in horizontal direction, ranging from -1.0 to 1.0")
    parser.add_argument('-mc', '--method_curv', type=str, default="disc",
                        help="Method for computing curv. Available methods: disc " +
                             "| vmtkfactor | vmtkit | spline")
    parser.add_argument('-ma', '--method_angle', type=str, default="plane",
                        help="Method for computing siphon angle. Available methods: plane | " +
                             "itplane | itplane_clip | maxcurv | smooth | discrete | frac | odrline | MISR ")
    parser.add_argument('-bd', '--boundary', nargs='+', default=None,
                        help='Boundary of grid, as a list: [alpha_min, alpha_max, beta_min, beta_max]')

    parser.add_argument("-g", "--grid-size", type=int, default=50,
                        help="Size of n x n matrix used for computing a set of discrete values")

    args = parser.parse_args()

    return dict(input_filepath=args.ifile, kappa=args.curvature, theta=args.angle,
                alpha=args.alpha, beta=args.beta, boundary=args.boundary, grid_size=args.grid_size,
                method_curv=args.method_curv, method_angle=args.method_angle, proj=args.proj)


if __name__ == "__main__":
    initialize(**read_command_line())

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

import operator
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from scipy import interpolate
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelextrema

# Local import
from morphman.common import *


def estimate_alpha_and_beta(input_filepath, quantity_to_compute, boundary, radius, grid_size, value_change,
                            method_angle,
                            method_curv, region_of_interest, region_points):
    """
    Imports a matrix of parameter values corresponding
    to a (alpha,beta) point and perform spline interpolation
    to make a parameter surface.
    Parameter surface is used to find intersection with
    initial parameter value plus / minus one standard deviation.
    Three criteria to find suggested alpha-beta value from intersection.

    Args:
        input_filepath (str): Surface model filename and location where data is stored.
        quantity_to_compute(str): Parameter name.
        boundary (list): Boundary of searching grid.
        radius (float): Minimum radius of circle to search outside of.
        grid_size (int): Size of searching grid ( grid_size x grid_size matrix)
        value_change (float): Desired change in curvature / bend angle to achieve
        method_angle (str): Method for computing angle.
        method_curv (str): Method for computing curvature.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'landmarking']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
    """
    # Get grid values
    base_path = get_path_names(input_filepath)

    # Get region points

    # Compute discrete values for quantity
    if type(boundary[-1]) is str:
        boundary = np.asarray(boundary, dtype=float)
    data = compute_quantities(input_filepath, boundary, quantity_to_compute, method_curv, method_angle,
                              region_of_interest, region_points, n=grid_size, projection=False)
    # Get grid boundary
    amin, amax, bmin, bmax = float(boundary[0]), float(boundary[1]), float(boundary[2]), float(boundary[3])

    # Set standard deviations used to find intersection

    # Defined SD planes for curvature
    # Tolerance added for adjusting SD
    # if there are no intersections found
    def value_plus(tolerance=0.0):
        return initial_value + value_change - tolerance

    def value_minus(tolerance=0.0):
        return initial_value - value_change + tolerance

    n = len(data)
    alpha = np.linspace(amin, amax, n)
    beta = np.linspace(bmin, bmax, n)
    alpha_long = np.linspace(amin, amax, 300)
    beta_long = np.linspace(bmin, bmax, 300)
    xx, yy = np.meshgrid(alpha, beta)

    points = np.zeros((n, 2))
    for i in range(len(xx)):
        points[i] = [alpha[i], beta[i]]

    # Spline interpolation
    f = interpolate.SmoothBivariateSpline(xx.ravel(), yy.ravel(), data.ravel())

    initial_value = f(0, 0)
    methods = [value_plus, value_minus]

    # Find intersecting points
    # Reduces SD if no points are found
    for plane in methods:
        zeros = alpha_beta_intersection(plane, f, alpha_long, beta_long)
        if len(zeros) == 0:
            empty = True

            # Leeway tolerance for matching quantity on interpolated surface
            tol = 0.005 if quantity_to_compute == "curvature" else 0.1
            max_iter = 50
            iterations = 0

            print("-- Found no points..Adjusting SD")
            while empty and iterations < max_iter:
                print("-- Iterations: %i" % (iterations + 1))
                zeros = alpha_beta_intersection(plane, f, alpha_long, beta_long, tol)
                if len(zeros) > 0:
                    empty = False
                iterations += 1
                if quantity_to_compute == "curvature":
                    tol += 0.001
                elif quantity_to_compute == "angle":
                    tol += 0.2

        # Check points and apply criteria
        # to find suggested values for alpha and beta
        if len(zeros) > 0:
            points = []
            for p in zeros:
                if (plane.__name__ == "value_plus" and quantity_to_compute == "curvature") \
                        or (plane.__name__ == "value_minus" and quantity_to_compute == "angle"):
                    if p[1] < 0:
                        points.append(p)
                else:
                    if p[1] > 0:
                        points.append(p)

            suggested_point = points[0]
            dist = 1e9
            for p in points[1:]:
                dist_tmp = la.norm(np.array(p))
                if radius < dist_tmp < dist:
                    dist = dist_tmp
                    suggested_point = p

            # Write points to file
            write_alpha_beta_point(base_path, suggested_point, plane.__name__, quantity_to_compute)


def compute_quantities(input_filepath, boundary, quantity, method_curv, method_angle, region_of_interest, region_points,
                       n=50, projection=False):
    """
    Initialization for computing curvature and angle.
    Values are either printed to terminal or stored in a (n x n) matrix.

    Args:
        input_filepath (str): Base location of surface model files.
        boundary (list): Bounds of grid for computing quantities.
        method_curv (str): Method used to compute curvature.
        method_angle (str): Method used to compute angle.
        quantity (string): Quantity to compute, either curvature or angle.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'landmarking']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        n (int): Determines matrix size when computing multiple values.
        projection (bool): Projects angle into 2d plane if True.

    Returns:
        values (ndarray): (n x n) matrix with quantity values
    """

    #  Initialize storage data
    values = np.zeros((n, n))
    amin, amax, bmin, bmax = boundary[0], boundary[1], boundary[2], boundary[3]
    alphas = np.linspace(amin, amax, n)
    betas = np.linspace(bmin, bmax, n)
    k = 0
    for i, alpha in enumerate(alphas):
        for j, beta in enumerate(betas):
            print("Iteration %i of %i" % (k + 1, n * n))
            if quantity == "curvature":
                value, _ = compute_curvature(input_filepath, alpha, beta, method_curv, None, False, region_of_interest,
                                             region_points)

            elif quantity == "angle":
                value, _ = compute_angle(input_filepath, alpha, beta, method_angle, None,
                                         region_of_interest, region_points, projection)
            values[i, j] = value
            k += 1

    return values


def compute_angle(input_filepath, alpha, beta, method, new_centerlines,
                  region_of_interest, region_points, projection=False):
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
        new_centerlines (vtkPolyData): New centerline.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'landmarking']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        projection (bool): True / False for computing 2D / 3D angle.

    Returns:
        new_deg (float): New angle of a vessel bend from a manipulated centerline.
    Returns:
        deg (float): Old angle of a vessel bend from a manipulated centerline.
    """
    # Get base path
    base_path = get_path_names(input_filepath)

    # Centerline path
    centerline_path = base_path + "_centerline.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Extract old centerline
    if not path.exists(centerline_path):

        # Compute centerlines
        inlet, outlets = get_centers(surface, base_path)

        print("-- Compute centerlines and Voronoi diagram")
        centerlines, _, _ = compute_centerlines(inlet, outlets, centerline_path,
                                                capped_surface, resampling=0.1,
                                                smooth=False, base_path=base_path)
    else:
        centerlines = read_polydata(centerline_path)

    # Get region of interest
    point_path = base_path + "_anterior_bend.particles"
    if region_of_interest == "landmarking":
        if not path.exists(point_path):
            raise RuntimeError(("The given .particles file: %s does not exist. Please run" +
                                " landmarking with automated_landmarking.py first.") % point_path)
        region_points = np.loadtxt(point_path)
    else:
        _, _, _, region_points,_ = get_line_to_change(capped_surface, centerlines,
                                                    region_of_interest, "bend", region_points, 0)
        region_points = [[region_points[3 * i], region_points[3 * i + 1], region_points[3 * i + 2]]
                         for i in range(len(region_points) // 3)]
    p1 = region_points[0]
    p2 = region_points[1]

    if new_centerlines is None:
        centerlines, new_centerlines = get_new_centerlines(centerlines, region_points, alpha, beta, p1, p2)

    # Get new siphon and prepare
    id1, id2, moved_id1, moved_id2, moved_p1, moved_p2 = get_moved_siphon(new_centerlines, centerlines, p1, p2)

    # Extract region of interest
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

    # Get direction to the point furthest away (dx)
    direction = "vertical"
    clipping_points_vtk = vtk.vtkPoints()
    for point in [p1, p2]:
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

        frenet_t1 = get_array("FrenetTangent", siphon, k=3)
        frenet_t2 = get_array("FrenetTangent", moved_siphon, k=3)

        p1_1, p1_id = find_closest_point(frenet_t1[-1], 0, max_id, p2, siphon)
        p2_2, p2_id = find_closest_point(frenet_t1[0], max_id, siphon.GetNumberOfPoints(), p1, siphon)

        newp1, np1_id = find_closest_point(frenet_t2[-1], 0, newmax_id, moved_p2, moved_siphon)
        newp2, np2_id = find_closest_point(frenet_t2[0], newmax_id,
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

    elif method == "smooth":
        allmaxcurv = argrelextrema(cutcurv, np.greater)[0]
        allnewmaxcurv = argrelextrema(newcutcurv, np.greater)[0]

        tmpcurv = cutcurv
        while len(allmaxcurv) > 2:
            tmpcurv = gaussian_filter(tmpcurv, 2)
            allmaxcurv = argrelextrema(tmpcurv, np.greater)[0]

        tmpnewcurv = newcutcurv
        while len(allnewmaxcurv) > 2:
            tmpnewcurv = gaussian_filter(tmpnewcurv, 2)
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

            deg = find_angle_odr(d1, d2, projection)
            new_deg = find_angle_odr(newd1, newd2, projection)

    elif method == "MISR":
        multiplier = 1.5
        n1 = siphon.GetNumberOfPoints()
        n2 = moved_siphon.GetNumberOfPoints()
        rad1 = siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(0)
        rad2 = siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(n1 - 1)
        newrad1 = moved_siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(0)
        newrad2 = moved_siphon.GetPointData().GetArray(radiusArrayName).GetTuple1(n2 - 1)

        pa, ra = move_past_sphere(siphon, p1, rad1, 0, step=1, stop=n1 - 1, X=multiplier)
        pb, rb = move_past_sphere(siphon, p2, rad2, n1 - 1, step=-1, stop=0, X=multiplier)
        new_pa, ra = move_past_sphere(moved_siphon, moved_p1, newrad1, 0, step=1, stop=n2 - 1, X=multiplier)
        new_pb, rb = move_past_sphere(moved_siphon, moved_p2, newrad2, n2 - 1, step=-1, stop=0, X=multiplier)

        deg, l1, l2 = find_angle(pa, pb, p1, p2, projection)
        new_deg, nl1, nl2 = find_angle(new_pa, new_pb, moved_p1, moved_p2, projection)

    else:
        if method == "frac":
            n_values = [5]
            left = [2]
            right = [3]
            i = 0
            dx = 1. / n_values[i]
            ida = int(id1 + (id2 - id1) * left[i] * dx)
            idb = int(id1 + (id2 - id1) * right[i] * dx)
            pa = siphon_splined.GetPoints().GetPoint(ida)
            pb = siphon_splined.GetPoints().GetPoint(idb)

            ida = int(moved_id1 + (moved_id2 - moved_id1) * left[i] * dx)
            idb = int(moved_id1 + (moved_id2 - moved_id1) * right[i] * dx)
            new_pa = moved_siphon.GetPoints().GetPoint(ida)
            new_pb = moved_siphon.GetPoints().GetPoint(idb)

            deg, l1, l2 = find_angle(pa, pb, p1, p2, projection)
            new_deg, nl1, nl2 = find_angle(new_pa, new_pb, moved_p1, moved_p2, projection)

        elif method in ["plane", "itplane", "itplane_clip", "maxcurv", "smooth",
                        "discrete", "maxdist"]:
            frac = 4. / 5.
            if method == "itplane_clip":
                id_mid = (p2_id - p1_id) / 2.
                new_id_mid = (np2_id - np1_id) / 2.
                if max_id > id_mid:
                    ida = int((max_id - p1_id) * frac)
                    idb = int((max_id - p1_id) * (1 + (1 - frac)))
                    pa = siphon.GetPoints().GetPoint(ida + p1_id)
                    pb = siphon.GetPoints().GetPoint(idb + p1_id)
                else:
                    idb = int((p2_id - max_id) * (1 + (1 - frac)))
                    ida = int((p2_id - max_id) * frac)
                    pa = siphon.GetPoints().GetPoint(ida)
                    pb = siphon.GetPoints().GetPoint(idb)

                if newmax_id > new_id_mid:
                    ida = int((newmax_id - np1_id) * frac)
                    idb = int((newmax_id - np1_id) * (1 + (1 - frac)))
                    new_pa = moved_siphon.GetPoints().GetPoint(ida + np1_id)
                    new_pb = moved_siphon.GetPoints().GetPoint(idb + np1_id)
                else:
                    ida = int((np2_id - newmax_id) * frac)
                    idb = int((np2_id - newmax_id) * (1 + (1 - frac)))
                    new_pa = moved_siphon.GetPoints().GetPoint(ida)
                    new_pb = moved_siphon.GetPoints().GetPoint(idb)

                deg, l1, l2 = find_angle(pa, pb, p1, p2, projection)
                new_deg, nl1, nl2 = find_angle(new_pa, new_pb, moved_p1, moved_p2, projection)

            else:
                ida = int(max_id * frac)
                idb = int(max_id * (1 + (1 - frac)))
                pa = siphon.GetPoints().GetPoint(ida)
                pb = siphon.GetPoints().GetPoint(idb)

                ida = int(newmax_id * frac)
                idb = int(newmax_id * (1 + (1 - frac)))
                new_pa = moved_siphon.GetPoints().GetPoint(ida)
                new_pb = moved_siphon.GetPoints().GetPoint(idb)

                deg, l1, l2 = find_angle(pa, pb, p1, p2, projection)
                new_deg, nl1, nl2 = find_angle(new_pa, new_pb, moved_p1, moved_p2, projection)

    return new_deg, deg


def compute_curvature(input_filepath, alpha, beta, method, new_centerlines, compute_original, region_of_interest,
                      region_points):
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
        compute_original (bool): Computes old curvature value if True.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'landmarking']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint

    Returns:
        new_maxcurv (float): Maximum curvature within the manipulated region of interest.
    Returns:
        old_maxcurv (float): Maximum curvature within the original region of interest.
    """
    # Input filename
    base_path = get_path_names(input_filepath)

    # Centerline filename
    centerline_path = base_path + "_centerline.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Extract old centerline
    if not path.exists(centerline_path):
        # Compute centerlines
        inlet, outlets = get_centers(surface, base_path)

        print("-- Compute centerlines and Voronoi diagram")
        centerlines, _, _ = compute_centerlines(inlet, outlets, centerline_path, capped_surface, resampling=0.1,
                                                smooth=False, base_path=base_path)
    else:
        centerlines = read_polydata(centerline_path)

    # Get region of interest
    point_path = base_path + "_anterior_bend.particles"
    if region_of_interest == "landmarking":
        if not path.exists(point_path):
            raise RuntimeError(("The given .particles file: %s does not exist. Please run" +
                                " landmarking with automated_landmarking.py first.") % point_path)
        region_points = np.loadtxt(point_path)
    else:
        _, _, _, region_points,_ = get_line_to_change(capped_surface, centerlines,
                                                    region_of_interest, "bend", region_points, 0)
        region_points = [[region_points[3 * i], region_points[3 * i + 1], region_points[3 * i + 2]]
                         for i in range(len(region_points) // 3)]
    p1 = region_points[0]
    p2 = region_points[1]

    if new_centerlines is None:
        print("-- Maniuplating centerline manually")
        centerlines, new_centerlines = get_new_centerlines(centerlines, region_points, alpha, beta, p1, p2)

    # Extract centerline points and ids
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
        new_curvature = gaussian_filter(curv_fac, 5)

        if compute_original:
            line_fac = vmtk_centerline_geometry(centerline, smooth=True, iterations=100, factor=factor)
            curv_fac = get_array("Curvature", line_fac)
            curvature = gaussian_filter(curv_fac, 5)

    # 2) VMTK - Iteration variance
    elif method == "vmtkit":
        it = 150
        line_it = vmtk_centerline_geometry(new_centerline, smooth=True, iterations=it, factor=1.0)
        curv_it = get_array("Curvature", line_it)
        new_curvature = gaussian_filter(curv_it, 5)

        if compute_original:
            line_it = vmtk_centerline_geometry(centerline, smooth=True, iterations=it, factor=1.0)
            curv_it = get_array("Curvature", line_it)
            curvature = gaussian_filter(curv_it, 5)

    # 3) Splines
    elif method == "spline":
        nknots = 50
        siphon_splined, siphon_curv = spline_centerline(new_centerline, get_curv=True,
                                                        isline=True, nknots=nknots)
        new_curvature = gaussian_filter(siphon_curv, 5)

        if compute_original:
            siphon_splined, siphon_curv = spline_centerline(centerline, get_curv=True,
                                                            isline=True, nknots=nknots)
            curvature = gaussian_filter(siphon_curv, 5)

    # 4) Default: Discrete derivatives
    elif method == "disc":
        neigh = 20
        line_di, curv_di = discrete_geometry(new_centerline, neigh=neigh)
        new_curvature = gaussian_filter(curv_di, 5)

        if compute_original:
            line_di, curv_di = discrete_geometry(centerline, neigh=neigh)
            curvature = gaussian_filter(curv_di, 5)

    old_maxcurv = max(curvature[id1 + 10:id2 - 10]) if compute_original else None
    new_maxcurv = max(new_curvature[id1_new + 10:id2_new - 10])

    return new_maxcurv, old_maxcurv


def get_new_centerlines(centerlines, region_points, alpha, beta, p1, p2):
    """
    Perform manual centerline manipulation in order to use the
    centerline as a proxy for computing a selected quantity (angle, curvature).
    Takes the original centerline, the points defining the region of interest and
    the manipulation factors alpha and beta before performing a manual movement of
    the centerline. Returns the new, manipulated centerline as a VTK object.

    Args:
        centerlines (vtkPolyData): Centerlines including diverging centerlines
        region_points (ndarray): List of region points
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        p1: First region point
        p2: Second region point

    Returns:
        centerlines (vtkPolyData): Centerlines excluding diverging centerlines
    Returns:
        new_centerlines (vtkPolyData): New centerlines including diverging centerlines
    """
    centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        find_region_of_interest_and_diverging_centerlines(centerlines, region_points)
    new_centerlines = centerlines
    diverging_id = None if len(diverging_ids) == 0 else diverging_ids[0]
    if beta != 0.0:
        direction = "horizont"
        middle_points, middle_ids = get_spline_points(extract_single_line(new_centerlines, 0), beta, direction,
                                                      region_points_vtk)
        dx_p1 = middle_points[0] - p1
        new_centerlines = move_centerlines(new_centerlines, dx_p1, p1, p2, diverging_id,
                                           diverging_centerlines, direction, merge_lines=True)
    if alpha != 0.0:
        direction = "vertical"
        middle_points, middle_ids, dx = get_spline_points(extract_single_line(new_centerlines, 0), alpha, direction,
                                                          region_points_vtk)
        merge_lines = False if beta != 0.0 else True
        new_centerlines = move_centerlines(new_centerlines, dx, p1, p2, diverging_id, diverging_centerlines, direction,
                                           merge_lines=merge_lines)

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
        id1 (int): ID of first clipping point.
    Returns:
        id2 (int): ID of second clipping point.
    Returns:
        moved_id1 (int): New ID of first clipping point.
    Returns:
        moved_id1 (int): New ID of secon clipping point.
    Returns:
        moved_p1 (ndarray): New position of first clipping point.
    Returns:
        moved_p2 (ndarray): New position of second ipping point.
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


def find_angle(pa, pb, p1, p2, projection):
    """
    Compute the angle between two vectors
    a = pA - p1 and b = pB - p2
    using the classical formula.

    Args:
        pa (ndarray): Point along the centerline.
        pb (ndarray): Point along the centerline.
        p1 (ndarray): Point along the centerline.
        p2 (ndarray): Point along the centerline.
        projection (bool): True / False for 2D / 3D angle.

    Returns:
        deg (float): Angle between vectors.
        vector_a (ndarray): First vector.
        vector_b (ndarraty): Second vector.
    """
    if not projection:
        vector_a = np.array([pa[0] - p1[0], pa[1] - p1[1], pa[2] - p1[2]])
        vector_b = np.array([pb[0] - p2[0], pb[1] - p2[1], pb[2] - p2[2]])
    else:
        vector_a = np.array([0, pa[1] - p1[1], pa[2] - p1[2]])
        vector_b = np.array([0, pb[1] - p2[1], pb[2] - p2[2]])
    costheta = (vector_a.dot(vector_b)) / (la.norm(vector_a) * la.norm(vector_b))
    angle = np.arccos(costheta)
    deg = (angle * 180 / np.pi)

    return deg, vector_a, vector_b


def find_angle_odr(d1, d2, projection):
    """
    Compute the angle between two vectors
    d1 and d2 using the classical formula.
    Used for the ODR-method, specifically.

    Args:
        d1 (ndarray): First vector
        d2 (ndarray): Second vector
        projection (bool): True / False for 2D / 3D angle.

    Returns:
        deg (float): Angle between vectors.
    """
    if d1.dot(d2) > 0:
        d1 = -d1
    if projection:
        d1[0] = 0
        d2[0] = 0

    costheta = (d1.dot(-d2)) / (la.norm(d1) * la.norm(-d2))
    angle = np.arccos(costheta)
    deg = (angle * 180 / np.pi)

    return deg, d1, d2


def write_alpha_beta_point(base_path, suggested_point, method, quantity_to_compute):
    """
    Write suggested choice of (alpha, beta) to file.

    Args:
        base_path (str): Path to file directory.
        suggested_point (ndarray): Array containing alpha and beta value
        method (str): Info about parameter, and increase / decrease of parameter.
        quantity_to_compute (str): Quantity to compute.
    """
    save_path = base_path + "_alphabeta_values.txt"
    alpha = suggested_point[0]
    beta = suggested_point[1]
    name = quantity_to_compute + "_" + method.split("_")[-1]
    with open(save_path, "a") as f:
        f.write("%s alpha=%.4f beta=%.4f \n" % (name, alpha, beta))


def alpha_beta_intersection(method, f, alphas, betas, tol=0.0):
    """
    Iterate through values of alpha and beta
    and find intersection points with a given tolerance
    between initial value and initial value plus/minus SD.

    Args:
        method (function): Plane defining initial value plus/minus SD.
        f (interp2d): Interpolated surface of curvature or angle values.
        alphas (ndarray): List of alpha values.
        betas (ndarray): List of beta values.
        tol (float): Tolerance used for adjusting SD if needed.

    Returns:
        zeros (ndarray): Array containing (alpha,beta) tuples which intersect.
    """
    zeros = []
    for i in alphas:
        for j in betas:
            diff = abs(f(j, i) - method(tol))
            if "c" in method.__name__:
                if diff < 0.001:
                    zeros.append([i, j])
            else:
                if diff < 0.05:
                    zeros.append([i, j])
    return zeros


def read_command_line():
    """
    Read arguments from commandline
    """
    description = "Algorithm used to compute the recommended value for " + \
                  "alpha and beta based on surface interpolation and a " + \
                  "given limit, depending on the quantity to be computed. " + \
                  "Primarily implemented for computing angle and curvature values." + \
                  "Computes selected quantities using the centerline as a proxy."

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required named arguments')

    # Required arguments
    required.add_argument('-i', '--ifile', type=str, default=None,
                          help="Path to the surface model", required=True)
    required.add_argument("-q", "--quantity", type=str, default="curvature",
                          help="Parameter to compute. Choose between 'curvature' and 'angle'", required=True,
                          choices=['curvature', 'angle'])
    # Optional arguments
    parser.add_argument("-r", "--region-of-interest", type=str, default="manual",
                        choices=["manual", "commandline", "landmarking"],
                        help="The method for defining the region to be changed. There are" +
                             " three options: 'manual', 'commandline', 'landmarking'. In" +
                             " 'manual' the user will be provided with a visualization of the" +
                             " input surface, and asked to provide an end and start point of the" +
                             " region of interest. Note that not all algorithms are robust over" +
                             " bifurcations. If 'commandline' is provided, then '--region-points'" +
                             " is expected to be provided. Finally, if 'landmarking' is" +
                             " given, it will look for the output from running" +
                             " automated_landmarking.py.")
    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" +
                             " argument have to be given. The method expects two points" +
                             " which defines the start and end of the region of interest. " +
                             " Example providing the points (1, 5, -1) and (2, -4, 3):" +
                             " --region-points 1 5 -1 2 -4 3")
    parser.add_argument('-rad', '--radius', type=float, default=0.15,
                        help="Radius of bounding circle, limiting the choice of alpha and beta")
    parser.add_argument('-b', '--boundary', nargs='+', default=[-0.2, 1, -0.2, 1],
                        help='Bounds of grid, as a list: [alpha_min, alpha_max, beta_min, beta_max]')
    parser.add_argument("-g", "--grid-size", type=int, default=25,
                        help="Size of n x n matrix used for computing a set of discrete values")
    parser.add_argument('-c', '--value-change', type=float, default=0.05,
                        help='Desired change in curvature / bend angle to achieve. Algorithm computes' +
                             'recommended values of alpha and beta for both plus and minus this change.')
    parser.add_argument('-mc', '--method_curv', type=str, default="disc",
                        help="Method for computing curv. Available methods: disc " +
                             "| vmtkfactor | vmtkit | spline",
                        choices=['disc', 'vmtkfactor', 'vmtkit', 'spline'])
    parser.add_argument('-ma', '--method_angle', type=str, default="plane",
                        help="Method for computing siphon angle. Available methods: plane | " +
                             "itplane | itplane_clip | maxcurv | smooth | discrete | frac | odrline | MISR ",
                        choices=['plane', 'itplane', 'itplane_clip', 'maxcurv', 'smooth', 'discrete', 'frac',
                                 'odrline', 'misr'])

    args = parser.parse_args()

    return dict(input_filepath=args.ifile, quantity_to_compute=args.quantity,
                radius=args.radius, boundary=args.boundary, grid_size=args.grid_size, value_change=args.value_change,
                method_angle=args.method_angle, method_curv=args.method_curv, region_points=args.region_points,
                region_of_interest=args.region_of_interest)


if __name__ == "__main__":
    estimate_alpha_and_beta(**read_command_line())

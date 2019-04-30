##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from scipy.interpolate import splrep, splev
from scipy.signal import resample

from morphman.common.common import *
from morphman.common.vessel_reconstruction_tools import create_parent_artery_patches
from morphman.common.vmtk_wrapper import *
from morphman.common.vmtkpointselector import vmtkPickPointSeedSelector
from morphman.common.vtk_wrapper import *


def get_bifurcating_and_diverging_point_data(centerline, centerline_bif, tol):
    """
    Locate bifurcating point and diverging points
    in a bifurcation.
    End points are set based on the MISR at
    the selected points.

    Args:
        centerline (vtkPolyData): Centerline from inlet to relevant outlets.
        centerline_bif (vtkPolyData): Centerline through bifurcation.
        tol (float): Tolerance parameter.

    Returns:
        data (dict): Contains info about diverging point locations.
    """
    # Sort centerline to start at inlet
    cl1 = extract_single_line(centerline, 0)
    cl2 = extract_single_line(centerline, 1)

    # Declear dictionary to hold results
    data = {"bif": {}, 0: {}, 1: {}}

    # Find lower clipping point
    n_points = min(cl1.GetNumberOfPoints(), cl2.GetNumberOfPoints())
    for i in range(0, n_points):
        point_0 = cl1.GetPoint(i)
        point_1 = cl2.GetPoint(i)
        distance_between_points = get_distance(point_0, point_1)
        if distance_between_points > tol:
            center = cl1.GetPoint(i)
            r = cl1.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
            break

    end, r_end, id_end = move_past_sphere(cl1, center, r, i, step=-1)
    data["bif"]["end_point"] = end
    data["bif"]["div_point"] = center

    # Find the diverging points for the bifurcation
    # continue further downstream in each direction and stop when
    # a point is closer than tol, then move point MISR * X
    locator = vtk_point_locator(centerline_bif)

    for counter, cl in enumerate([cl1, cl2]):
        for i in range(i, cl.GetNumberOfPoints(), 1):
            tmp_point = cl.GetPoint(i)
            closest_point_id = locator.FindClosestPoint(tmp_point)
            closest_point = centerline_bif.GetPoint(closest_point_id)
            distance_between_points = get_distance(tmp_point, closest_point)
            if distance_between_points < tol * 4:
                center = cl.GetPoint(i)
                r = cl1.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
                break

        end, r_end, id_end = move_past_sphere(cl, center, r, i, step=1,
                                              stop=i * 100, scale_factor=1)
        data[counter]["end_point"] = end
        data[counter]["div_point"] = center

    return data


def get_manipulated_centerlines(patch_cl, dx, p1, p2, diverging_id, diverging_centerlines, direction, merge_lines=True):
    """Given a centerline (patch_cl), move the centerline a distance (dx) between two
    points (p1 and p2).

    Args:
        patch_cl (vtkPolyData): Centerlines excluding diverging centerlines.
        dx (ndarray): Direction to move geometry.
        p1 (vtkPolyData): First region point.
        p2: (vtkPolyData): Second region point.
        diverging_id (int): List of index where centerlines diverge from region of interest.
        diverging_centerlines (vtkPolyData): Centerlines which diverge from region of interest.
        direction (str): Manipulation direction parameter.
        merge_lines (bool): Merge centerlines and diverging centerlines.

    Returns:
        centerline (vtkPolyData): Manipulated centerline.
    """
    if diverging_id is not None and merge_lines:
        patch_cl = vtk_merge_polydata([patch_cl, diverging_centerlines])

    number_of_points = patch_cl.GetNumberOfPoints()
    number_of_cells = patch_cl.GetNumberOfCells()

    centerline = vtk.vtkPolyData()
    centerline_points = vtk.vtkPoints()
    centerline_cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, number_of_points)

    count = 0
    for i in range(number_of_cells):
        line = extract_single_line(patch_cl, i)
        centerline_cell_array.InsertNextCell(line.GetNumberOfPoints())

        radius_array_data = line.GetPointData().GetArray(radiusArrayName).GetTuple1

        locator = vtk_point_locator(line)
        id1 = locator.FindClosestPoint(p1)
        if diverging_id is not None and i == (number_of_cells - 1):
            # Note: Reuse id2 and idmid from previous loop
            pass
        else:
            id2 = locator.FindClosestPoint(p2)
            idmid = int((id1 + id2) * 0.5)

        for p in range(line.GetNumberOfPoints()):
            point = line.GetPoint(p)
            cl_id = locator.FindClosestPoint(point)

            if direction == "horizont":
                if diverging_id is not None and i == (number_of_cells - 1) and diverging_id < cl_id:
                    dist = -dx * (diverging_id - idmid) ** 0.5 / (diverging_id - idmid) ** 0.5
                else:
                    if cl_id < id1:
                        dist = dx
                    elif id1 <= cl_id < idmid:
                        dist = dx * (idmid ** 2 - cl_id ** 2) / (idmid ** 2 - id1 ** 2)
                    elif idmid <= cl_id < (id2 - 1):
                        dist = -dx * (cl_id - idmid) ** 0.5 / (id2 - idmid) ** 0.5
                    else:
                        dist = -dx

            elif direction == "vertical":
                if diverging_id is not None and i == (number_of_cells - 1) and diverging_id < cl_id:
                    dist = 4 * dx * (diverging_id - id1) * (id2 - diverging_id) / (id2 - id1) ** 2
                else:
                    if id1 <= cl_id <= id2:
                        dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2
                    else:
                        dist = 0

            point = np.asarray(point)
            centerline_points.InsertNextPoint(point + dist)
            radius_array.SetTuple1(count, radius_array_data(p))
            centerline_cell_array.InsertCellPoint(count)
            count += 1

    centerline.SetPoints(centerline_points)
    centerline.SetLines(centerline_cell_array)
    centerline.GetPointData().AddArray(radius_array)

    return centerline


def get_centerline_between_clipping_points(centerline_relevant_outlets, data):
    """Get the centerline between two clipping points.

    Args:
        centerline_relevant_outlets (vtkPolyData): Centerline to the two relevant outlets.
        data (dict): A dictionary data.

    Returns:
        centerline (vtkPolyData): Return clipped centerline.
    """
    line0 = extract_single_line(centerline_relevant_outlets, 0)
    line1 = extract_single_line(centerline_relevant_outlets, 1)
    lines = []
    for i, line in enumerate([line0, line1]):
        loc = vtk_point_locator(line)
        tmp_id_dau = loc.FindClosestPoint(data[i]["end_point"])
        tmp_id_bif = loc.FindClosestPoint(data["bif"]["end_point"])
        lines.append(extract_single_line(line, 0, start_id=tmp_id_bif, end_id=tmp_id_dau))

    centerline = vtk_merge_polydata(lines)

    return centerline


def get_diverging_point_id(centerline1, centerline2, tol):
    """
    Find ID of diverging point;
    where two input centerlines diverge
    due to a bifurcation.

    Args:
        centerline1 (vtkPolyData): First centerline.
        centerline2 (vtkPolyData): Second centerline.
        tol (float): Tolerance.

    Returns:
        i (int): ID at diverging point.
    """
    # Find clipping points
    n_points = min(centerline1.GetNumberOfPoints(), centerline2.GetNumberOfPoints())
    get_point1 = centerline1.GetPoints().GetPoint
    get_point2 = centerline2.GetPoints().GetPoint

    for i in range(0, n_points):
        distance_between_points = get_distance(get_point1(i), get_point2(i))
        if distance_between_points > tol:
            break

    return i


def get_curvilinear_coordinate(line):
    """
    Get curvilinear coordinates along
    an input centerline.

    Args:
        line (vtkPolyData): Input centerline

    Returns:
        curv_coor (ndarray): Array of abscissa points.
    """
    curv_coor = np.zeros(line.GetNumberOfPoints())
    for i in range(line.GetNumberOfPoints() - 1):
        pnt1 = np.asarray(line.GetPoints().GetPoint(i))
        pnt2 = np.asarray(line.GetPoints().GetPoint(i + 1))
        curv_coor[i + 1] = np.sqrt(np.sum((pnt1 - pnt2) ** 2)) + curv_coor[i]

    return curv_coor


def get_centerline_tolerance(centerline, n=50):
    """
    Finds tolerance based on
    average length between first N points
    along the input centerline.

    Args:
        centerline (vtkPolyData): Centerline data.
        n (int): Number of points

    Returns:
        tolerance (float): Tolerance value.
    """
    line = extract_single_line(centerline, 0)
    length = get_curvilinear_coordinate(line)
    tolerance = np.mean(length[1:n] - length[:n - 1]) / 2.0

    return tolerance


def get_clipped_diverging_centerline(centerline, clip_start_point, clip_end_id):
    """
    Clip the opthamlic artery if present.

    Args:
        centerline (vtkPolyData): Line representing the opthalmic artery centerline.
        clip_start_point (tuple): Point at entrance of opthalmic artery.
        clip_end_id (int): ID of point at end of opthalmic artery.

    Returns:
        patch_eye (vtkPolyData): Voronoi diagram representing opthalmic artery.
    """
    points = [clip_start_point, centerline.GetPoint(clip_end_id)]
    div_points = vtk.vtkPoints()
    for p in points:
        div_points.InsertNextPoint(p)

    patch_cl = create_parent_artery_patches(centerline, div_points, siphon=True)

    return patch_cl


def get_line_to_change(surface, centerline, region_of_interest, method, region_points,
                       stenosis_length):
    """
    Extract and spline part of centerline
    within the geometry where
    area variations will be performed.

    Args:
        surface (vtkPolyData): Surface model.
        centerline (vtkPolyData): Centerline in geometry.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'first_line']
        method (str): Determines which kind of manipulation is performed.
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint.
        stenosis_length (float): Multiplier used to determine the length of the stenosis-affected area.

    Returns:
        line_to_change (vtkPolyData): Part of centerline.
    """
    if region_of_interest == "first_line":
        tol = get_centerline_tolerance(centerline)
        line2 = extract_single_line(centerline, 0)
        number_of_points2 = line2.GetNumberOfPoints()

        # Iterate through lines and find diverging point
        n = centerline.GetNumberOfLines()
        point_ids = []
        for j in range(1, n):
            line1 = extract_single_line(centerline, j)
            number_of_points1 = line1.GetNumberOfPoints()

            n_points = min(number_of_points1, number_of_points2)
            for i in range(n_points):
                point1 = line1.GetPoints().GetPoint(i)
                point2 = line2.GetPoints().GetPoint(i)
                if get_distance(point1, point2) > tol:
                    point_id = i
                    break
            point_ids.append(point_id)

        start_id = 0
        end_id = min(point_ids)
        cl_id = point_ids.index(end_id)

        region_points = list(line2.GetPoint(start_id)) + list(line2.GetPoint(end_id))

    elif region_of_interest in ["commandline", "landmarking", "manual"]:
        # Get points from the user
        if region_of_interest == "manual":
            print("\nPlease select region of interest in the render window.")
            stenosis_point_id = vtk.vtkIdList()
            first = True
            while stenosis_point_id.GetNumberOfIds() not in [1, 2]:
                if not first:
                    print("Please provide only one or two points, try again")

                # Select point on surface
                seed_selector = vmtkPickPointSeedSelector()
                seed_selector.SetSurface(surface)
                if method == "variation" or method == "area":
                    seed_selector.text = "Press space to select the start and endpoint of the" + \
                                         " region of interest, 'u' to undo.\n"
                elif method == "stenosis":
                    seed_selector.text = "Press space to select, the center of a new" + \
                                         " stenosis (one point),\nOR place two points on each side" + \
                                         " of an existing stenosis to remove it, \'u\' to undo."
                elif method == "bend":
                    seed_selector.text = "Press space to select the start and end of the" + \
                                         " bend that you want to manipulate, 'u' to undo.\n"
                seed_selector.Execute()
                stenosis_point_id = seed_selector.GetTargetSeedIds()
                first = False

            region_points = []
            for i in range(stenosis_point_id.GetNumberOfIds()):
                region_points += surface.GetPoint(stenosis_point_id.GetId(i))

            print("-- The chosen region points:", region_points)

        # Get locator
        locator = vtk_point_locator(centerline)

        if len(region_points) == 3:
            # Project point onto centerline
            region_points = centerline.GetPoint(locator.FindClosestPoint(region_points))
            point1 = region_points

            # Get relevant line
            tol = get_centerline_tolerance(centerline)
            cl_id = -1
            dist = 1e10
            while dist > tol / 10:
                cl_id += 1
                line = extract_single_line(centerline, cl_id)
                tmp_loc = vtk_point_locator(line)
                tmp_id = tmp_loc.FindClosestPoint(point1)
                dist = get_distance(point1, line.GetPoint(tmp_id))

            # Get length of stenosis
            misr = get_point_data_array(radiusArrayName, line)
            length = stenosis_length * misr[tmp_loc.FindClosestPoint(point1)]

            # Get ids of start and stop
            centerline_length = get_curvilinear_coordinate(line)
            center = centerline_length[tmp_id]
            region_of_interest_id = (center - length <= centerline_length) * (centerline_length <= center + length)
            start_id = np.argmax(region_of_interest_id)
            end_id = region_of_interest_id.shape[0] - 1 - np.argmax(region_of_interest_id[::-1])

        else:
            point1 = region_points[:3]
            point2 = region_points[3:]
            point1 = centerline.GetPoint(locator.FindClosestPoint(point1))
            point2 = centerline.GetPoint(locator.FindClosestPoint(point2))
            region_points[:3] = point1
            region_points[3:] = point2

            distance1 = []
            distance2 = []
            ids1 = []
            ids2 = []
            for i in range(centerline.GetNumberOfLines()):
                line = extract_single_line(centerline, i)
                tmp_loc = vtk_point_locator(line)
                ids1.append(tmp_loc.FindClosestPoint(point1))
                ids2.append(tmp_loc.FindClosestPoint(point2))
                cl_point1 = line.GetPoint(ids1[-1])
                cl_point2 = line.GetPoint(ids2[-1])
                distance1.append(get_distance(point1, cl_point1))
                distance2.append(get_distance(point2, cl_point2))

            tol = get_centerline_tolerance(centerline) / 10
            total_distance = (np.array(distance1) < tol) * (np.array(distance2) < tol)
            cl_id = np.argmax(total_distance)

            if total_distance[cl_id] == 0:
                raise RuntimeError("The two points provided have to be on the same " +
                                   " line (from inlet to outlet), and not at two different" +
                                   " outlets")
            start_id = min(ids1[cl_id], ids2[cl_id])
            end_id = max(ids1[cl_id], ids2[cl_id])

    # Extract and spline a single line
    line_to_change = extract_single_line(centerline, cl_id, start_id=start_id, end_id=end_id)

    remaining_centerlines = []
    diverging_centerlines = []
    start_point = line_to_change.GetPoint(0)
    end_point = line_to_change.GetPoint(line_to_change.GetNumberOfPoints() - 1)
    tol = get_centerline_tolerance(centerline) * 4
    for i in range(centerline.GetNumberOfLines()):
        line = extract_single_line(centerline, i)
        locator = vtk_point_locator(line)
        id1 = locator.FindClosestPoint(start_point)
        id2 = locator.FindClosestPoint(end_point)
        p1_tmp = line.GetPoint(id1)
        p2_tmp = line.GetPoint(id2)
        close_start = get_distance(start_point, p1_tmp) < tol
        close_end = get_distance(end_point, p2_tmp) < tol

        # Check if the centerline is going through both or None of the start and end of
        # the region of interest
        if close_start == close_end:
            if close_start:
                if start_id != 0:
                    tmp = extract_single_line(centerline, i, start_id=0, end_id=start_id - 1)
                    remaining_centerlines.append(tmp)
                tmp = extract_single_line(centerline, i, start_id=end_id + 1,
                                          end_id=line.GetNumberOfPoints() - 1)
                remaining_centerlines.append(tmp)
            else:
                remaining_centerlines.append(line)
        else:
            diverging_centerlines.append(line)

    # Find diverging ids
    diverging_ids = []
    main_line = extract_single_line(centerline, 0)
    id_start = 0
    for line in diverging_centerlines:
        id_end = min([line.GetNumberOfPoints(), main_line.GetNumberOfPoints()])
        for i in np.arange(id_start, id_end):
            p_div = np.asarray(line.GetPoint(i))
            p_cl = np.asarray(main_line.GetPoint(i))
            if get_distance(p_div, p_cl) > tol * 10:
                diverging_ids.append(i)
                break

    # Spline the single line
    nknots = min(line_to_change.GetNumberOfPoints() // 2, 25)
    line_to_change = compute_splined_centerline(line_to_change, nknots=nknots, isline=True)

    if len(diverging_centerlines) == 0:
        diverging_centerlines = None
    remaining_centerlines = vtk_merge_polydata(remaining_centerlines)

    return line_to_change, remaining_centerlines, diverging_centerlines, region_points, diverging_ids


def get_region_of_interest_and_diverging_centerlines(centerlines_complete, region_points):
    """Extract the centerline between the region points.

    Args:
        centerlines_complete (vktPolyData): Complete set of centerlines in geometry.
        region_points (ndarray): Two points determining the region of interest.

    Returns:
        centerlines (vtkPolyData): Centerlines excluding divering lines.
        diverging_centerlines (vtkPolyData): Centerlines diverging from the region of interest.
        region_points (ndarray): Sorted region points.
        region_points_vtk (vtkPoints): Sorted region points as vtkData.
        diverging_ids (list): List of indices where a diverging centerline starts.
    """
    centerlines = []
    diverging_centerlines = []
    p1 = region_points[0]
    p2 = region_points[1]

    # Search for divering centerlines
    tol = get_centerline_tolerance(centerlines_complete) * 4
    for i in range(centerlines_complete.GetNumberOfLines()):
        line = extract_single_line(centerlines_complete, i)
        locator = vtk_point_locator(line)
        id1 = locator.FindClosestPoint(p1)
        id2 = locator.FindClosestPoint(p2)
        p1_tmp = line.GetPoint(id1)
        p2_tmp = line.GetPoint(id2)
        if get_distance(p1, p1_tmp) < tol and get_distance(p2, p2_tmp) < tol:
            centerlines.append(line)
        else:
            diverging_centerlines.append(line)

    # Sort and set clipping points to vtk object
    centerline = centerlines[0]
    locator = vtk_point_locator(centerline)
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    if id1 > id2:
        region_points = region_points[::-1]
        id1, id2 = id2, id1

    region_points_vtk = vtk.vtkPoints()
    for point in np.asarray(region_points):
        region_points_vtk.InsertNextPoint(point)

    # Find diverging point(s)
    diverging_ids = []
    for line in diverging_centerlines:
        id_end = min([line.GetNumberOfPoints(), centerline.GetNumberOfPoints()])
        for i in np.arange(id1, id_end):
            p_div = np.asarray(line.GetPoint(i))
            p_cl = np.asarray(centerline.GetPoint(i))
            if get_distance(p_div, p_cl) > tol:
                diverging_ids.append(i)
                break

    centerlines = vtk_merge_polydata(centerlines)
    diverging_centerlines = vtk_merge_polydata(diverging_centerlines) if len(diverging_centerlines) > 0 else None

    return centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids


def compute_discrete_derivatives(line, neigh=10):
    """Compute the curvature and torsion of a line using 'neigh' number of neighboring points.

    Args:
        line (vtkPolyData): Line to compute geometry from.
        neigh (int): Number of naboring points.

    Returns:
        line (vtkPolyData): Output line with geometrical parameters.
        curv (vtkPolyData): Output line with geometrical parameters.
    """
    n_points = line.GetNumberOfPoints()

    # Compute cumulative chord length
    t = np.zeros(n_points)
    p = []
    for i in range(n_points):
        p.append(np.array(list(line.GetPoint(i))))
        p[i] = np.array(p[i])

    norms = [la.norm(p[j] - p[j - 1]) for j in range(1, n_points)]
    s = sum(norms)
    for i in range(1, n_points):
        s1 = sum(norms[:i + 1])
        t[i] = s1 / s

    # Radius of sliding neighbourhood
    m = neigh

    dxdt = np.zeros(n_points)
    dydt = np.zeros(n_points)
    dzdt = np.zeros(n_points)

    x = np.zeros(n_points)
    y = np.zeros(n_points)
    z = np.zeros(n_points)

    for i in range(n_points):
        x[i] = p[i][0]
        y[i] = p[i][1]
        z[i] = p[i][2]

    for i in range(0, m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(0, 2 * m + 1)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(0, 2 * m + 1)]) / t_sum

    for i in range(m, n_points - m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(i - m, i + m + 1)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(i - m, i + m + 1)]) / t_sum

    for i in range(n_points - m, n_points):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(n_points - 2 * m, n_points)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum

    dgammadt = []
    dgammadt_norm = np.zeros(n_points)
    for i in range(n_points):
        dgammadt.append(np.array([dxdt[i], dydt[i], dzdt[i]]))
        dgammadt_norm[i] = la.norm(dgammadt[i])

    tg = []
    for i in range(n_points):
        tg.append(dgammadt[i] / dgammadt_norm[i])

    t1 = np.zeros(n_points)
    t2 = np.zeros(n_points)
    t3 = np.zeros(n_points)

    for i in range(n_points):
        t1[i] = tg[i][0]
        t2[i] = tg[i][1]
        t3[i] = tg[i][2]

    dt1dt = np.zeros(n_points)
    dt2dt = np.zeros(n_points)
    dt3dt = np.zeros(n_points)

    for i in range(0, m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(0, 2 * m + 1)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(0, 2 * m + 1)]) / t_sum

    for i in range(m, n_points - m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(i - m, i + m + 1)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(i - m, i + m + 1)]) / t_sum

    for i in range(n_points - m, n_points):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(n_points - 2 * m, n_points)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum

    dtgdt = []
    dtgdt_norm = np.zeros(n_points)
    for i in range(n_points):
        dtgdt.append(np.array([dt1dt[i], dt2dt[i], dt3dt[i]]))
        dtgdt_norm[i] = la.norm(dtgdt[i])

    curv = np.zeros(n_points)
    for i in range(n_points):
        curv[i] = dtgdt_norm[i] / dgammadt_norm[i]

    curv = resample(curv, n_points)

    return line, curv


def get_k1k2_basis(curvature, line):
    """
    Create a k1-k2 basis used to determine
    the location of each bend of the carotid
    siphon.

    Args:
        curvature (floats): Curvature array.
        line (vtk): Centerline points.

    Returns:
        line (vtk): Centerline points including k1-k2 basis.
    """

    # The ParallelTransportNormals and the FrenetTangent is not orthonormal
    # (but close) from vmtk. Using the GramSchmidt proses gives E2 and fixes
    # the non-orthogonality
    E1 = get_point_data_array("ParallelTransportNormals", line, k=3)
    T = get_point_data_array("FrenetTangent", line, k=3)
    E2 = np.zeros((E1.shape[0], 3))

    for i in range(E1.shape[0]):
        V = np.eye(3)
        V[:, 0] = T[i, :]
        V[:, 1] = E1[i, :]
        V = gram_schmidt(V)

        E1[i, :] = V[:, 1]
        E2[i, :] = V[:, 2]

    # Compute k_1, k_2 furfilling T' = curv(s)*N(s) = k_1(s)*E_1(s) + k_2(s)*E_2(s).
    # This is simply a change of basis for the curvature vector N. The room
    # of k_1 and k_2 can be used to express both the curvature and the
    # torsion.
    N = get_point_data_array("FrenetNormal", line, k=3)

    k2 = (curvature.T * (E1[:, 1] * N[:, 0] - N[:, 1] * E1[:, 0]) /
          (E2[:, 1] * E1[:, 0] - E2[:, 0] * E1[:, 1]))[0]
    k1 = (-(curvature.T * N[:, 0] + k2 * E2[:, 0]) / E1[:, 0])[0]

    for k in [(k1, "k1"), (k2, "k2")]:
        k_array = create_vtk_array(k[0], k[1])
        line.GetPointData().AddArray(k_array)

    return line


def compute_splined_centerline(line, get_curv=False, isline=False, nknots=50, get_stats=True, get_misr=True):
    """
    Given the knots and coefficients of a B-spline representation,
    evaluate the value of the smoothing polynomial and its derivatives.
    This is a wrapper around the FORTRAN routines splev and splder of FITPACK.

    Args:
        line (vtkPolyData): Centerline points.
        get_curv (bool): Computes curvature profile if True.
        isline (bool): Determines if centerline object is a line or points.
        nknots (int): Number of knots.
        get_stats (bool): Determines if curve attribuites are computed or not.
        get_misr (bool): Determines if MISR values are computed or not.

    Returns:
        line (vtkPolyData): Splined centerline data.
    Returns:
        curv (ndarray): Curvature profile.
    """

    if not isline:
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
    if get_misr:
        misr = get_point_data_array(radiusArrayName, line)

    curv_coor = get_curvilinear_coordinate(line)
    for i in range(data.shape[0]):
        data[i, :] = line.GetPoint(i)

    t = np.linspace(curv_coor[0], curv_coor[-1], nknots + 2)[1:-1]
    fx = splrep(curv_coor, data[:, 0], k=4, t=t)
    fy = splrep(curv_coor, data[:, 1], k=4, t=t)
    fz = splrep(curv_coor, data[:, 2], k=4, t=t)

    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    if get_misr:
        data = np.zeros((len(curv_coor), 4))
        data[:, 3] = misr[:, 0]
        header = ["X", "Y", "Z", radiusArrayName]
    else:
        data = np.zeros((len(curv_coor), 3))
        header = ["X", "Y", "Z"]

    data[:, 0] = fx_
    data[:, 1] = fy_
    data[:, 2] = fz_

    line = convert_numpy_data_to_polydata(data, header)

    # Let vmtk compute curve attributes
    if get_stats:
        line = vmtk_compute_geometric_features(line, smooth=False)

    if get_curv:
        # Analytical curvature
        dlsfx = splev(curv_coor, fx, der=1)
        dlsfy = splev(curv_coor, fy, der=1)
        dlsfz = splev(curv_coor, fz, der=1)

        ddlsfx = splev(curv_coor, fx, der=2)
        ddlsfy = splev(curv_coor, fy, der=2)
        ddlsfz = splev(curv_coor, fz, der=2)

        c1xc2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
        c1xc2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
        c1xc2_3 = ddlsfy * dlsfx - ddlsfx * dlsfy

        curvature = np.sqrt(c1xc2_1 ** 2 + c1xc2_2 ** 2 + c1xc2_3 ** 2) / (dlsfx ** 2 + dlsfy ** 2 + dlsfz ** 2) ** 1.5

        return line, curvature
    else:
        return line

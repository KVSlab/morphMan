def centerline_div(centerline1, centerline2, tol):
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
        distance_between_points = distance(get_point1(i), get_point2(i))
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


def get_tolerance(centerline, n=50):
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
    tolerance = np.mean(length[1:n] - length[:n - 1]) / divergingRatioToSpacingTolerance

    return tolerance

def get_data(centerline, centerline_bif, tol):
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
        distance_between_points = distance(point_0, point_1)
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
    locator = get_locator(centerline_bif)

    for counter, cl in enumerate([cl1, cl2]):
        for i in range(i, cl.GetNumberOfPoints(), 1):
            tmp_point = cl.GetPoint(i)
            closest_point_id = locator.FindClosestPoint(tmp_point)
            closest_point = centerline_bif.GetPoint(closest_point_id)
            distance_between_points = distance(tmp_point, closest_point)
            if distance_between_points < tol * 4:
                center = cl.GetPoint(i)
                r = cl1.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
                break

        end, r_end, id_end = move_past_sphere(cl, center, r, i, step=1,
                                              stop=i * 100, X=1)
        data[counter]["end_point"] = end
        data[counter]["div_point"] = center

    return data


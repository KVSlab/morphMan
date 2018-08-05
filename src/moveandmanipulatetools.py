import operator
import numpy.linalg as la
from common import *
from argparse import ArgumentParser
from scipy.signal import argrelextrema, resample
from scipy.ndimage.filters import gaussian_filter as gauss

# Local import
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *


def get_clipping_points(dirpath, filename):
    """
    Read clipping points from file

    Args:
        dirpath (str): Location of directory.
        filename (str): Name of clipping point file.

    Returns:
        clipping_points (ndarray): Clipping points.
    """
    particles = path.join(dirpath, filename)
    all_points = np.loadtxt(particles)
    clipping_points = all_points

    return clipping_points


def connect_line(line):
    """
    Create edges between points defining a
    centerline, used to construct a
    discrete line in 3D.

    Args:
        line (vtkPoints): Points defining a centerline.

    Returns:
        line (vtkPolyData): Discrete centerline data.
    """
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


def move_line_horizontally(patch_cl, ID1, ID2, dx_p1, clip=False, eye=False, side=None):
    """
    Iterate through centerline and move line based on a profile
    for horizontal movement. Includes special treatment of
    opthalmic artery if present.

    Args:
        patch_cl (vtkPolyData): Centerline data.
        ID1 (int): Index of first clipping point.
        ID2 (int): Index of second clipping point.
        dx_p1 (ndarray): Direction to move upstream.
        clip (bool): Determines which part of geometry is being moved, True if siphon.
        eye (bool): Determines presence of opthamlic artery.
        side (str): Determines location relative to the middle of the siphon.

    Returns:
        newline (vtkPolyData): Manipulated centerline.
    """

    centerline_loc = get_locator(patch_cl)
    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if clip:
        if eye:
            l1 = extract_single_line(patch_cl, 0)
            l2 = extract_single_line(patch_cl, 1)
            l3 = extract_single_line(patch_cl, 2)
            test_cl = merge_data([l1, l2, l3])

            ID1 = 0
            ID2 = len(get_curvilinear_coordinate(test_cl))
            idmid = int((ID1 + ID2) / 2.)
        else:
            ID1 = 0
            ID2 = len(get_curvilinear_coordinate(patch_cl))
            idmid = int((ID1 + ID2) / 2.)

        for p in range(patch_cl.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(patch_cl.GetPoint(p))

            if cl_id < idmid:
                dist = dx_p1 * (idmid ** 2 - cl_id ** 2) / (idmid ** 2 - ID1 ** 2)
            else:
                if cl_id <= (ID2 - 1):
                    dist = -dx_p1 * (cl_id - idmid) ** (0.5) / (ID2 - idmid) ** (0.5)
                else:
                    locator = get_locator(test_cl)
                    pp = patch_cl.GetPoint(cl_id)
                    id_main = locator.FindClosestPoint(pp)
                    dist = -dx_p1 * (id_main - idmid) ** (0.5) / (id_main - idmid) ** (0.5)

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


def move_points_vertically(line, dx):
    """
    Iterate through centerline points and move line based on a profile
    for vertical movement.

    Args:
        line (vtkPolyData): Centerline data.
        dx (ndarray): Direction to move centerline.

    Returns:
        newline (vtkPolyData): Manipulated centerline.
    """

    centerline_loc = get_locator(line)

    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    ID1 = 0
    ID2 = len(get_curvilinear_coordinate(line))
    for p in range(line.GetNumberOfPoints()):
        cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))

        dist = 4 * dx * (cl_id - ID1) * (ID2 - cl_id) / (ID2 - ID1) ** 2

        points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
        verts.InsertNextCell(1)
        verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(line.GetPointData().GetArray(radiusArrayName))
    return newline


def move_line_vertically(line, dx, ID1_0, clip_ID=None, eye=False):
    """
    Iterate through centerline and move line based on a profile
    for vertical movement. Includes special treatment of
    opthalmic artery if present.

    Args:
        line (vtkPolyData): Centerline data.
        dx (ndarray): Direction to move vertically.
        ID1_0 (int): Index of first clipping point.
        clip_ID (int): Index of opthalmic artery entrance if present.
        eye (bool): Determines presence of opthamlic artery.

    Returns:
        newline (vtkPolyData): Manipulated centerline.
    """
    centerline_loc = get_locator(line)

    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if eye:
        l1 = extract_single_line(line, 0)
        l2 = extract_single_line(line, 1)
        l3 = extract_single_line(line, 2)
        test_cl = merge_data([l1, l2, l3])

        ID1 = I1 = 0
        ID2 = len(get_curvilinear_coordinate(test_cl))
        IDmid = int((ID1 + ID2) / 2.)
        I2 = ID2 - 1

        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))

            if cl_id <= I2:
                dist = 4 * dx * (cl_id - I1) * (I2 - cl_id) / (I2 - I1) ** 2
            else:
                cl_id = clip_ID - ID1_0 + int((ID2 - (clip_ID - ID1_0)) * 0.4)
                dist = 4 * dx * (cl_id - ID1) * (ID2 - cl_id) / (ID2 - ID1) ** 2

            points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)
    else:
        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(line))

        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))
            dist = 4 * dx * (cl_id - ID1) * (ID2 - cl_id) / (ID2 - ID1) ** 2

            points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(line.GetPointData().GetArray(radiusArrayName))

    return newline


def move_perp(n, P, Z, alpha):
    """
    Find directions for manipulation
    in the vertical direction.

    Args:
        n (ndarray): Normal vector to plane through clipping points.
        P (ndarray): Clipping points.
        Z (ndarray): Points along the centerline.
        alpha (float): Extension / Compression factor.

    Returns:
        dZ (ndarray): Directions to points along the  centerline.
        dx (ndarray): Direction to move the centerline.
    """

    p1 = P[0]
    p2 = P[1]

    # Find midpoint and point furthest away
    dist = []
    for z in Z:
        d = la.norm(np.cross((z - p1), (z - p2))) / la.norm(p2 - p1)
        dist.append(d)

    D_id = dist.index(max(dist))
    D = max(dist)
    z_m = Z[D_id]

    # Vector from line to Z_max and projection onto plane
    v = (z_m - p2) - (z_m - p2).dot(p2 - p1) * (p2 - p1) / la.norm(p2 - p1) ** 2
    PV = v - v.dot(n) * n

    # Find distances
    P1 = (z_m - p1).dot(p2 - p1) * (p2 - p1) / la.norm(p2 - p1) ** 2
    P1 = p1 + P1
    V = P1 + v
    PV1 = P1 + PV

    # Move points
    dZ = []
    for i in range(len(Z)):
        dz = np.array(PV) * dist[i] / D * alpha
        dZ.append(Z[i] + dz)

    dx = (PV1 - P1) * alpha

    return dZ, dx


def move_para(n, P, Z, beta):
    """
    Find directions for manipulation
    in the horizontal direction.

    Args:
        n (ndarray): Normal vector to plane through clipping points.
        P (ndarray): Clipping points.
        Z (ndarray): Points along the centerline.
        beta (float): Extension / Compression factor.

    Returns:
        dZ (ndarray): Directions to points along the  centerline.
        zp_min (ndarray): Translation direction in upstream direction.
        zm_min (ndarray): Translation direction in downstream direction.
    """
    p1 = P[0]
    p2 = P[1]

    # Find normal from q
    q = [(p1 + p2) / 2.]
    qp2 = np.array(p2 - q)[0]
    qp1 = np.array(p1 - q)[0]
    s = q[0] - np.cross(qp2, n) * 3

    # Split points based on orientation
    # to q normal
    Z_p = []
    Z_p_dist = []
    Z_m = []
    Z_m_dist = []
    for z in Z:
        d = la.norm(np.cross((z - s), (z - q[0]))) / la.norm(q[0] - s)
        c = np.cross(s - q, z - q)
        if c[0][0] >= 0:
            Z_p.append(z)
            Z_p_dist.append(d)
        else:
            Z_m.append(z)
            Z_m_dist.append(d)

    # Move points
    D = la.norm(p1 - p2) / 2.
    dZ = []
    dZ.append(p1 + qp1 * beta)
    for i in range(len(Z_p)):
        dz = qp1 * Z_p_dist[i] / D * beta
        dZ.append(Z_p[i] + dz)

    for i in range(len(Z_m)):
        dz = qp2 * Z_m_dist[i] / D * beta
        dZ.append(Z_m[i] + dz)
    dZ.append(p2 + qp2 * beta)

    # Check if moved in right direction
    d_0 = la.norm(np.cross(Z_p[0] - q, Z_p[0] - s)) / la.norm(s - q)
    d_1 = la.norm(np.cross(dZ[1] - q, dZ[1] - s)) / la.norm(s - q)
    if d_1 < d_0:
        # Split points based on orientation
        # to q normal
        Z_p = []
        Z_p_dist = []
        Z_m = []
        Z_m_dist = []
        for z in Z:
            d = la.norm(np.cross((z - s), (z - q[0]))) / la.norm(q[0] - s)
            c = -np.cross(s - q, z - q)
            if c[0][0] >= 0:
                Z_p.append(z)
                Z_p_dist.append(d)
            else:
                Z_m.append(z)
                Z_m_dist.append(d)

        # Move points
        D = la.norm(p1 - p2) / 2.
        dZ = []
        dZ.append(p1 + qp1 * beta)
        for i in range(len(Z_p)):
            dz = qp1 * Z_p_dist[i] / D * beta
            dZ.append(Z_p[i] + dz)

        for i in range(len(Z_m)):
            dz = qp2 * Z_m_dist[i] / D * beta
            dZ.append(Z_m[i] + dz)
        dZ.append(p2 + qp2 * beta)

    zpid = Z_p_dist.index(min(Z_p_dist))
    zp_min = Z_p[zpid]
    zmid = Z_m_dist.index(min(Z_m_dist))
    zm_min = Z_m[zmid]

    return dZ, zp_min, zm_min


def best_plane(Z, P):
    """
    Find the least squares plane through
    the points in P and approximating the points
    in Z.

    Args:
        Z (ndarray): Array of points to approximate.
        P (ndarray): Array of points used as constraints.

    Returns:
        n (ndarray): Normal vector to plane.
    """
    # Defined matrices
    b = np.ones(len(Z))
    d = np.ones(len(P))

    # Create complete matrix
    ATA = np.transpose(Z).dot(Z)
    M0 = np.c_[ATA, np.transpose(P)]
    M1 = np.c_[P, np.zeros((len(P), len(P)))]
    M = np.r_[M0, M1]
    Y = np.r_[np.transpose(Z).dot(b), d]

    # Solve system
    x = la.solve(M, Y)
    a = x[0]
    b = x[1]
    c = x[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    # Define plane
    xmin = min(Z, key=operator.itemgetter(1))[0] - 4
    xmax = max(Z, key=operator.itemgetter(1))[0] + 4
    ymin = min(Z, key=operator.itemgetter(1))[1] - 4
    ymax = max(Z, key=operator.itemgetter(1))[1] + 4
    xx, yy = np.meshgrid(np.linspace(xmin, xmax, 15), np.linspace(ymin, ymax, 15))
    zz = (1 - a * xx - b * yy) / float(c)

    return n


def find_closest_point(dx, start, stop, P0, line):
    """
    Find point located closest to a given point P0.
    Searching from start to stop along the centerline.

    Args:
        dx (ndarray): Direction to search for point furthest away.
        start (int): Index to start searching.
        stop (int): Index to stop searching.
        P0 (ndarray): Point to search from.
        line (vtkPolyData): Centerline to search along.

    Returns:
        minP (ndarray): Point located closest to P0.
        minID (int): ID of point located closest to P0.
    """
    a = dx[0]
    b = dx[1]
    c = dx[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    # Define plane
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 100
    xx, yy = np.meshgrid(np.linspace(xmin, xmax, 150), np.linspace(ymin, ymax, 150))
    d = a * P0[0] + b * P0[1] + c * P0[2]
    zz = (d - a * xx - b * yy) / float(c)

    points = []
    for i in range(start, stop):
        p = line.GetPoint(i)
        points.append(np.array(p))

    dist_list = []
    for i, pcl in enumerate(points):
        v = pcl - np.array(P0)
        dist = abs(v.dot(n))
        dist_list.append(dist)

    minID = dist_list.index(min(dist_list)) + start
    minP = points[minID - start]

    return minP, minID


def find_furthest_points(dx, line):
    """
    Find point located furthes away from the line
    spanned of the clipping points p1 and p2.

    Args:
        dx (ndarray): Direction to search for point furthest away.
        line (vtkPolyData): Centerline to search along.

    Returns:
        maxP (ndarray): Point located furthest away.
        maxID (int): ID of point located furthest away.
    """
    P0 = line.GetPoint(0)
    a = dx[0]
    b = dx[1]
    c = dx[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    # Define plane
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 100
    xx, yy = np.meshgrid(np.linspace(xmin, xmax, 150), np.linspace(ymin, ymax, 150))
    d = a * P0[0] + b * P0[1] + c * P0[2]
    zz = (d - a * xx - b * yy) / float(c)

    points = []
    for i in range(line.GetNumberOfPoints()):
        p = line.GetPoint(i)
        points.append(np.array(p))

    dist_list = []
    for i, pcl in enumerate(points):
        v = pcl - np.array(P0)
        dist = abs(v.dot(n))
        dist_list.append(dist)

    maxID = dist_list.index(max(dist_list))
    maxP = points[maxID]

    return maxP, maxID


def get_spline_points(line, param, direction, clip_points):
    """
    Pick n uniformly selected points along the
    centerline from point P1 to P2, and move them.

    Args:
        line (vtkPolyData): Longest centerline in geometry.
        param (float): Extension / Compression factor.
        direction (str): Direction to move centerline.
        clip_points (vtkPoints): Clipping points.

    Returns:
        dz (ndarray): Points along the centerline.
        ids (ndarray): IDs of points along centerline.
        dx (ndarray): Direction to move geometry.
    """
    locator = get_locator(line)
    p1 = clip_points.GetPoint(0)
    p2 = clip_points.GetPoint(1)
    ID1 = locator.FindClosestPoint(p1)
    ID2 = locator.FindClosestPoint(p2)
    ID_mid = int((ID1 + ID2) / 2.)
    P = [p1, p2]

    # Select n uniformly spaced points
    n = 10
    points = []
    ids = np.zeros(n)
    dx = 1 / (n + 1.)
    for i in range(1, n + 1):
        ID = int(ID1 + (ID2 - ID1) * i * dx)
        ids[i - 1] = ID
        p = line.GetPoints().GetPoint(ID)
        points.append(np.array([p[0], p[1], p[2]]))

    for i in range(len(P)):
        P[i] = np.array([P[i][0], P[i][1], P[i][2]])

    n = best_plane(points, P)

    if direction == "vertical":
        dz, dx = move_perp(n, P, points, param)
        return dz, ids, dx

    elif direction == "horizont":
        dz, zp, zm = move_para(n, P, points, param)
        return dz, ids


def find_diverging_centerlines(centerlines, end_point):
    """
    Collect centerlines diverging from the longest
    centerline. 

    Args:
        centerlines (vtkPolyData): Collection of all centerlines.
        end_point (tuple): End point of relevant centerline.

    Returns:
        div_ids (list): ID where lines diverege.
    Returns:
        div_points (list): Points where lines diverge.
    Returns:
        centerlines (vtkPolyData): Collection of centerlines not diverging.
    Returns:
        div_lines (list): Collection of divering centerlines.  
    """
    # Start with longest line
    longest = extract_single_line(centerlines, 0)
    longest = vmtk_centerline_resampling(longest, 0.1)
    longest_locator = get_locator(longest)
    longest_end_id = longest_locator.FindClosestPoint(end_point)

    # Separate lines and divering lines
    lines = [longest]
    div_lines = []
    div_ids = []
    div_points = []
    n = centerlines.GetNumberOfCells() - 1
    tol = 0.40

    # Find diverging lines
    for i in range(n):
        div = False
        line_tmp = extract_single_line(centerlines, i + 1)
        line_tmp = vmtk_centerline_resampling(line_tmp, 0.1)
        stop_id = len(get_curvilinear_coordinate(line_tmp))
        if longest_end_id < stop_id:
            stop_id = longest_end_id
        for j in np.arange(stop_id):
            p_cl = np.asarray(longest.GetPoint(j))
            p_tmp = np.asarray(line_tmp.GetPoint(j))
            dist = distance(p_cl, p_tmp)
            if dist > tol and j < (longest_end_id-20):
                div_lines.append(line_tmp)
                div_ids.append(j)
                div_points.append(p_tmp)
                div = True
                break
        if not div:
            lines.append(line_tmp)

    centerlines = merge_data(lines)

    return div_ids, div_points, centerlines, div_lines


def clip_eyeline(eyeline, clip_start_point, clip_end_ID):
    """
    Clip the opthamlic artery if present.

    Args:
        eyeline (vtkPolyData): Line representing the opthalmic artery centerline.
        clip_star_point (tuple): Point at entrance of opthalmic artery.
        cip_end_ID (int): ID of point at end of opthalmic artery.

    Returns:
        patch_eye (vtkPolyData): Voronoi diagram representing opthalmic artery.
    """

    points = [clip_start_point, eyeline.GetPoint(clip_end_ID)]
    eye_points = vtk.vtkPoints()
    for p in points:
        eye_points.InsertNextPoint(p)

    patch_eye = CreateParentArteryPatches(eyeline, eye_points, siphon=True)

    return patch_eye


def find_ophthalmic_artery(centerlines, clip_pts):
    """
    Method checks if the geometry includes the opthamlic artery.
    Extracts the opthalmic artery if present, and determines its position.

    Args:
        centerlines (vtkPolyData): Complete set of centerlines.
        clip_pts (vtkPoints): Clipping points.

    Returns:
        eye (bool): True if opthamlic artery is present.
        clip_ID (long): ID where opthamlic artery is located along centerline.
        centerlines (vtkPolyData): Complete set of centerlines excluding opthamlic artery.
        eyeline (vtkPolyData): Centerline leading to opthalmic artery.
    """

    # Extract lines:
    lines = []
    n = centerlines.GetNumberOfCells()
    for i in range(n):
        lines.append(extract_single_line(centerlines, i))

    longest = lines[0]
    tol = 0.40

    # Find start and stop IDs along clipped curve
    locator_longest = get_locator(longest)
    p1 = clip_pts[0]
    p2 = clip_pts[1]
    ID1 = locator_longest.FindClosestPoint(p1)
    ID2 = locator_longest.FindClosestPoint(p2)

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
            p_eye = np.asarray(line.GetPoint(i))
            p_cl = np.asarray(longest.GetPoint(i))
            dist = la.norm(p_eye - p_cl)
            if dist > tol:
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


def get_vtk_clipping_points(line, clipping_points):
    """
    Store clipping points as VTK objects.
    Extract points as tuples and corresponding IDs.

    Args:
        line (vtkPolyData): Line representing longest single centerline.
        clipping_points (ndarray): Array containing two clpping points.

    Returns:
        p1 (tuple): First clipping point.
        p2 (tuple): Second clipping point.
        ID1 (long): ID of first clipping point.
        ID2 (long): ID of second clipping point.
        clip_points (vtkPoints): VTK objects containing the clipping points.
        clipping_points (ndarray): Array containing two clipping points.
    """
    locator = get_locator(line)
    ID1 = locator.FindClosestPoint(clipping_points[0])
    ID2 = locator.FindClosestPoint(clipping_points[1])
    if ID1 > ID2:
        clipping_points = clipping_points[::-1]
        ID1, ID2 = ID2, ID1

    # Set clipping points
    div_points = np.asarray(clipping_points)
    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)
    clip_points = points

    p1 = clip_points.GetPoint(0)
    p2 = clip_points.GetPoint(1)

    return p1, p2, ID1, ID2, clip_points, clipping_points

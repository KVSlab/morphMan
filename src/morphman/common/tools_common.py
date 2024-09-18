##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.
## This software is distributed WITHOUT ANY WARRANTY; without even
## the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE. See the above copyright notices for more information.

import json
from os import path

import numpy as np
import numpy.linalg as la
import vtk

from morphman.common.vtk_wrapper import get_vtk_array, get_vtk_point_locator


def get_path_names(input_filepath):
    """Takes the input folder path as argument, and returns the name of the case name, and
    the path to the parent directory

    Args:
        input_filepath (str): Input filepath

    Returns:
        base_path (str): Path to the surface, but without the extension
    """

    surface_name = input_filepath.split(path.sep)[-1].split(".")[0]
    folder = path.dirname(input_filepath)
    base_path = path.join(folder, surface_name)

    return base_path


def get_distance(point1, point2):
    """Distance between two points.

    Args:
        point1 (ndarray): A point
        point2 (ndarray): A point

    Returns:
        distance (float): Distance between point1 and point2
    """
    return np.sqrt(np.sum((np.asarray(point1) - np.asarray(point2)) ** 2))


def gram_schmidt(V):
    """Gram schmidt process of each column

    Args:
        V (numpy.array): A (n x n) matrix

    Returns:
        E (numpy.array): A (n x n) matrix where all the columns are orthogonal
    """
    V = 1.0 * V
    U = np.copy(V)

    def proj(u, v):
        return u * np.dot(v, u) / np.dot(u, u)

    for i in range(1, V.shape[1]):
        for j in range(i):
            U[:, i] -= proj(U[:, j], V[:, i])

    # normalize column
    den = (U**2).sum(axis=0) ** 0.5
    E = U / den

    return E


def get_parameters(folder):
    """Read the parameters in the info file.

    Args:
        folder (str): Path to folder.

    Returns:
        data (dict): The data in the info file.
    """
    # If info.json does not exist, return an empty dict
    if not path.isfile(folder + "_info.json"):
        return {}

    # Get dictionary
    with open(folder + "_info.json", "r") as infile:
        data = json.load(infile)

    return data


def write_parameters(data, folder):
    """Get the old parameters, then write the new parameters in data.

    Args:
        data (dict): New data to write to parameters
        folder (str): Path to data location.
    """
    # Get old parameters
    parameters = get_parameters(folder)

    # Add new parameters (can overwrite old as well)
    parameters.update(data)

    # Write data
    with open(folder + "_info.json", "w") as outfile:
        json.dump(parameters, outfile)


def convert_numpy_data_to_polydata(data, header, TNB=None, PT=None):
    """Converting a range of data to a vtk array.

    Args:
        data (numpy.ndarray): Data array.
        header (list): A list of names for each array.
        TNB (numpy.ndarray): Data array.
        PT (numpy.ndarray): Data array.

    Returns:
        line (vtkPolyData): Line-couple with all the new data.
    """
    line = vtk.vtkPolyData()
    cell_array = vtk.vtkCellArray()
    cell_array.InsertNextCell(data.shape[0])
    line_points = vtk.vtkPoints()

    info_array = []
    for i in range(3, data.shape[1]):
        radius_array = get_vtk_array(header[i], 1, data.shape[0])
        info_array.append(radius_array)

    if TNB is not None:
        for i in range(3):
            radius_array = get_vtk_array(header[i + data.shape[1]], 3, data.shape[0])
            info_array.append(radius_array)

    if PT is not None:
        start = data.shape[1] if TNB is None else data.shape[1] + 3
        for i in range(2):
            radius_array = get_vtk_array(header[i + start], 3, PT[0].shape[0])
            info_array.append(radius_array)

    for i in range(data.shape[0]):
        cell_array.InsertCellPoint(i)
        line_points.InsertNextPoint(data[i, :3])
        for j in range(3, data.shape[1]):
            info_array[j - 3].SetTuple1(i, data[i, j])

    if TNB is not None:
        for i in range(data.shape[0]):
            for j in range(data.shape[1] - 3, data.shape[1], 1):
                tnb_ = TNB[j - data.shape[1]][i, :]
                info_array[j].SetTuple3(i, tnb_[0], tnb_[1], tnb_[2])

    if PT is not None:
        start = data.shape[1] - 3 if TNB is None else data.shape[1]
        for i in range(PT[-1].shape[0]):
            for j in range(start, start + 2, 1):
                pt_ = PT[j - start][i, :]
                info_array[j].SetTuple3(i, pt_[0], pt_[1], pt_[2])

    line.SetPoints(line_points)
    line.SetLines(cell_array)
    for i in range(len(header) - 3):
        line.GetPointData().AddArray(info_array[i])

    return line


def get_sorted_outlets(outlets, outlet1, outlet2, dirpath):
    """
    Sort all outlets of the geometry given the two relevant outlets

    Args:
        outlets (list): List of outlet center points.
        outlet1 (list): Point representing first relevant oultet.
        outlet2 (list): Point representing second relevant oultet.
        dirpath (str): Location of info file.

    Returns:
        outlets (list): List of sorted outlet center points.
        outlet1 (list): Point representing first relevant oultet.
        outlet2 (list): Point representing second relevant oultet.
    """
    tmp_outlets = np.array(outlets).reshape(len(outlets) // 3, 3)
    outlet1_index = np.argsort(np.sum((tmp_outlets - outlet1) ** 2, axis=1))[0]
    outlet2_index = np.argsort(np.sum((tmp_outlets - outlet2) ** 2, axis=1))[0]
    tmp_outlets = tmp_outlets.tolist()

    if max(outlet1_index, outlet2_index) == outlet1_index:
        outlet1 = tmp_outlets.pop(outlet1_index)
        outlet2 = tmp_outlets.pop(outlet2_index)
    else:
        outlet2 = tmp_outlets.pop(outlet2_index)
        outlet1 = tmp_outlets.pop(outlet1_index)
    outlet_rest = (np.array(tmp_outlets).flatten()).tolist()
    outlets = outlet1 + outlet2 + outlet_rest
    data = {}

    for i in range(len(outlets) // 3):
        data["outlet" + str(i)] = outlets[3 * i : 3 * (i + 1)]

    write_parameters(data, dirpath)

    return outlets, outlet1, outlet2


def get_vertical_direction_parameters(n, region_points, cl_points, alpha):
    """
    Find directions for manipulation
    in the vertical direction.

    Args:
        n (ndarray): Normal vector to plane through clipping points.
        region_points (ndarray): Region points.
        cl_points (ndarray): Points along the centerline.
        alpha (float): Extension / Compression factor.

    Returns:
        dZ (ndarray): Directions to points along the  centerline.
        dx (ndarray): Direction to move the centerline.
    """
    # Find midpoint and point furthest away
    p1 = region_points[0]
    p2 = region_points[1]
    dist = []
    for point in cl_points:
        d = la.norm(np.cross((point - p1), (point - p2))) / la.norm(p2 - p1)
        dist.append(d)

    d_id = dist.index(max(dist))
    max_dist = max(dist)
    z_m = cl_points[d_id]

    # Vector from line to Z_max and projection onto plane
    v = (z_m - p2) - (z_m - p2).dot(p2 - p1) * (p2 - p1) / la.norm(p2 - p1) ** 2
    dv = v - v.dot(n) * n

    # Find distances
    dp = p1 + (z_m - p1).dot(p2 - p1) * (p2 - p1) / la.norm(p2 - p1) ** 2
    dpv = dp + dv

    # Move points
    z_plus_dz = []
    for i, _ in enumerate(cl_points):
        dz = np.array(dv) * dist[i] / max_dist * alpha
        z_plus_dz.append(cl_points[i] + dz)

    dx = (dpv - dp) * alpha

    return z_plus_dz, dx


def get_horizontal_direction_parameters(n, region_points, cl_points, beta):
    """
    Find directions for manipulation
    in the horizontal direction.

    Args:
        n (ndarray): Normal vector to plane through clipping points.
        region_points (ndarray): Clipping points.
        cl_points (ndarray): Points along the centerline.
        beta (float): Extension / Compression factor.

    Returns:
        dZ (ndarray): Directions to points along the  centerline.
        zp_min (ndarray): Translation direction in upstream direction.
    """
    p1 = region_points[0]
    p2 = region_points[1]

    # Find normal from q
    q = [(p1 + p2) / 2.0]
    qp2 = np.array(p2 - q)[0]
    qp1 = np.array(p1 - q)[0]
    s = q[0] - np.cross(qp2, n) * 3

    # Split points based on orientation to q normal
    points_p = []
    points_p_dist = []
    points_m = []
    points_m_dist = []
    for point in cl_points:
        d = la.norm(np.cross((point - s), (point - q[0]))) / la.norm(q[0] - s)
        c = np.cross(s - q, point - q)
        if c[0][0] >= 0:
            points_p.append(point)
            points_p_dist.append(d)
        else:
            points_m.append(point)
            points_m_dist.append(d)

    # Move points
    mid_point = la.norm(p1 - p2) / 2.0
    moved_points = [p1 + qp1 * beta]
    for i, _ in enumerate(points_p):
        dz = qp1 * points_p_dist[i] / mid_point * beta
        moved_points.append(points_p[i] + dz)

    for i, _ in enumerate(points_m):
        dz = qp2 * points_m_dist[i] / mid_point * beta
        moved_points.append(points_m[i] + dz)
    moved_points.append(p2 + qp2 * beta)

    # Check if moved in right direction
    d_0 = la.norm(np.cross(points_p[0] - q, points_p[0] - s)) / la.norm(s - q)
    d_1 = la.norm(np.cross(moved_points[1] - q, moved_points[1] - s)) / la.norm(s - q)
    if d_1 < d_0:
        # Split points based on orientation to q normal
        points_p = []
        points_p_dist = []
        points_m = []
        points_m_dist = []
        for z in cl_points:
            d = la.norm(np.cross((z - s), (z - q[0]))) / la.norm(q[0] - s)
            c = -np.cross(s - q, z - q)
            if c[0][0] >= 0:
                points_p.append(z)
                points_p_dist.append(d)
            else:
                points_m.append(z)
                points_m_dist.append(d)

        # Move points
        moved_points = [p1 + qp1 * beta]
        for i in range(len(points_p)):
            dz = qp1 * points_p_dist[i] / mid_point * beta
            moved_points.append(points_p[i] + dz)

        for i in range(len(points_m)):
            dz = qp2 * points_m_dist[i] / mid_point * beta
            moved_points.append(points_m[i] + dz)
        moved_points.append(p2 + qp2 * beta)

    zpid = points_p_dist.index(min(points_p_dist))
    zp_min = points_p[zpid]

    return moved_points, zp_min


def compute_least_square_plane(cl_points, region_points):
    """
    Find the least squares plane through
    the points in P and approximating the points
    in Z.

    Args:
        cl_points (ndarray): Array of points to approximate.
        region_points (ndarray): Array of points used as constraints.

    Returns:
        n (ndarray): Normal vector to plane.
    """
    # Defined matrices
    b = np.ones(len(cl_points))
    d = np.ones(len(region_points))

    # Create complete matrix
    ata = np.transpose(cl_points).dot(cl_points)
    m0 = np.c_[ata, np.transpose(region_points)]
    m1 = np.c_[region_points, np.zeros((len(region_points), len(region_points)))]
    m = np.r_[m0, m1]
    y = np.r_[np.transpose(cl_points).dot(b), d]

    # Solve system
    x = la.solve(m, y)
    a, b, c = x[0], x[1], x[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    return n


def get_closest_point(dx, start, stop, p0, line):
    """
    Find point located closest to a given point P0.
    Searching from start to stop along the centerline.

    Args:
        dx (ndarray): Direction to search for point closest.
        start (int): Index to start searching.
        stop (int): Index to stop searching.
        p0 (ndarray): Point to search from.
        line (vtkPolyData): Centerline to search along.

    Returns:
        minP (ndarray): Point located closest to P0.
        minID (int): ID of point located closest to P0.
    """
    a, b, c = dx[0], dx[1], dx[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    points = []
    for i in range(start, stop):
        p = line.GetPoint(i)
        points.append(np.array(p))

    dist_list = []
    for i, pcl in enumerate(points):
        v = pcl - np.array(p0)
        dist = abs(v.dot(n))
        dist_list.append(dist)

    min_id = dist_list.index(min(dist_list)) + start
    min_p = points[min_id - start]

    return min_p, min_id


def get_most_distant_point(dx, line):
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
    p0 = line.GetPoint(0)
    a, b, c = dx[0], dx[1], dx[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    points = []
    for i in range(line.GetNumberOfPoints()):
        p = line.GetPoint(i)
        points.append(np.array(p))

    dist_list = []
    for i, pcl in enumerate(points):
        v = pcl - np.array(p0)
        dist = abs(v.dot(n))
        dist_list.append(dist)

    max_id = dist_list.index(max(dist_list))
    max_p = points[max_id]

    return max_p, max_id


def get_direction_parameters(line, param, direction, clip_points):
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
    locator = get_vtk_point_locator(line)
    p1 = clip_points.GetPoint(0)
    p2 = clip_points.GetPoint(1)
    id1 = locator.FindClosestPoint(p1)
    id2 = locator.FindClosestPoint(p2)
    region_points = [p1, p2]

    for i, _ in enumerate(region_points):
        region_points[i] = np.array(
            [region_points[i][0], region_points[i][1], region_points[i][2]]
        )

    # Select n uniformly spaced points
    n = 10
    points = []
    ids = np.zeros(n)
    dx = 1 / (n + 1.0)
    for i in range(1, n + 1):
        id_ = int(id1 + (id2 - id1) * i * dx)
        ids[i - 1] = id_
        p = line.GetPoints().GetPoint(id_)
        points.append(np.array([p[0], p[1], p[2]]))

    n = compute_least_square_plane(points, region_points)

    if direction == "vertical":
        dz, dx = get_vertical_direction_parameters(n, region_points, points, param)
        return dz, ids, dx

    elif direction == "horizont":
        dz, zp = get_horizontal_direction_parameters(n, region_points, points, param)
        return dz, ids


def get_rotation_matrix(u, angle):
    """
    Get three dimensional rotation matrix based on Euler-Rodrigues formula

    Args:
        u (ndarray): Normal vector corresponding to rotation axis
        angle (float): Angle to rotate in radians

    Returns:
        rotation_matrix (ndarray): Rotation matrix
    """
    u_cross_matrix = np.asarray([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    u_outer = np.outer(u, u)
    rotation_matrix = (
        np.cos(angle) * np.eye(3)
        + np.sin(angle) * u_cross_matrix
        + (1 - np.cos(angle)) * u_outer
    )

    return rotation_matrix


def get_angle(a, b):
    """
    Compute angle between two lines,
    expressed as vectors.

    Args:
        a (ndarray): First vector
        b (ndarray): Second vector

    Returns:
        float: Angle between lines in radians
    """
    angle = np.arccos(np.dot(a, b) / (la.norm(a) * la.norm(b)))

    return angle

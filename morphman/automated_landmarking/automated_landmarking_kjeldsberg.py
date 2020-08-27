##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

import os

# Local import
from morphman.automated_landmarking.automated_landmarking_tools import *


def landmarking_kjeldsberg(centerline, base_path, smoothing_factor, iterations, smooth_line, resampling_step,
                           coronal_axis, mark_diverging_arteries_manually):
    """
    Perform automated classification of the internal carotid artery into seven segments
    C1 to C7. If model is too short, a minimum of four (C7 - C4) segments are located, adjusted accordingly.

    Interfaces between bends are mainly detected using the curvature profile of the centerline, and adjusted
    according to angles between normal vector planes at interfaces.

    Args:
        centerline (vtkPolyData): ICA centerline
        base_path (str): Path to current model
        smoothing_factor (float): Relaxation parameter in Laplacian smoothing
        iterations (int): Number of smoothing iterations
        smooth_line (boolean): Determines if ICA centerline is smoothed
        resampling_step (float): Resampling step
        coronal_axis (str) : Axis determining coronal coordinate
        mark_diverging_arteries_manually (boolean): Mark Ophthalmic and Posterior communicating artery manually
    """
    short_model = False
    very_short_model = False

    # Centerline computations
    centerlines_path = base_path + "_centerline.vtp"
    centerline_complete = read_polydata(centerlines_path)
    if resampling_step is not None:
        centerline = vmtk_resample_centerline(centerline, length=resampling_step)
        centerline_complete = vmtk_resample_centerline(centerline_complete, length=resampling_step)

    # Compute geometric properties from a smoothed representation of the centerline
    line = vmtk_compute_centerline_attributes(centerline)
    line = vmtk_compute_geometric_features(line, smooth_line, factor=smoothing_factor, iterations=iterations)
    curvature = get_point_data_array("Curvature", line)
    tangent = get_point_data_array("FrenetTangent", line, k=3)
    curvilinear_coordinates = get_curvilinear_coordinate(line)

    # Find diverging arteries if they exist
    if mark_diverging_arteries_manually:
        classify_with_diverging_arteries, ophthalmic_id, p_com_a_id = mark_diverging_arteries(centerline_complete, line)
    else:
        classify_with_diverging_arteries, ophthalmic_id, p_com_a_id = find_diverging_centerlines(centerline_complete,
                                                                                                 line)

    # Remove additional noise
    curvature = gaussian_filter(curvature, 15)

    # Set coronal coordinate and find C4-C5 interface
    coordinates = get_centerline_coordinates(line, curvilinear_coordinates)
    coronal_coordinates = coordinates[coronal_axis]

    c4_c5 = find_c4_c5_interface(curvilinear_coordinates, coronal_coordinates)

    # Find max and min curvature
    min_point_ids = list(argrelextrema(curvature, np.less)[0])
    max_point_ids = list(argrelextrema(curvature, np.greater)[0])

    # Find interface C3-C4
    c3_c4, start, very_short_model = find_c3_c4_interface(c4_c5, max_point_ids, min_point_ids, very_short_model)

    # Find interface C2-C3
    if not very_short_model:
        c2_c3, short_model, start = find_c2_c3_interface(max_point_ids, min_point_ids, start, short_model,
                                                         very_short_model)
    else:
        short_model = True

    if not short_model:
        # Check if C3-C4 and C2-C3 interfaces bends
        c2_c3 = adjust_c2_c3_interface(c2_c3, c3_c4, tangent)

        # Find interface C1-C2
        c1_c2 = find_c1_c2_interface(max_point_ids, min_point_ids, start)

        # Check if C1-C2 and C2-C3 interfaces bends
        c1_c2 = adjust_c1_c2_interface(c1_c2, c2_c3, tangent)

        # Move C2-C3 interfaces towards first downstream bend
        c2_c3_closest_max = max_point_ids[(max_point_ids > c2_c3).nonzero()[0][0]]
        c2_c3 = int((c2_c3 + c2_c3_closest_max) * 0.5)

        # Move C1-C2 interfaces towards first downstream bend
        c1_c2_closest_max = max_point_ids[(max_point_ids > c1_c2).nonzero()[0][0]]
        c1_c2 = int((c1_c2 + c1_c2_closest_max) * 0.5)

    # Set c5-c6 interface when C5 wedge is 30 degrees
    c5_c6 = find_c5_c6_interface(c4_c5, tangent)

    # Find interface C6-C7
    c6_c7 = find_c6_c7_interface(c5_c6, max_point_ids, curvilinear_coordinates)

    # Adjust C5-C6 interface: length of C5 < 2/3 * C6
    c5_c6 = adjust_c5_c6_interface(c4_c5, c5_c6, c6_c7, curvilinear_coordinates)

    # Adjust C6-C7 interface: length of C7 < C6
    c6_c7 = adjust_c6_c7_interface(c5_c6, c6_c7, curvilinear_coordinates)

    if classify_with_diverging_arteries:
        c5_length = c5_c6 - c4_c5

        # Overwrite c5-c6 and c6-c7 interface based on arteries
        c5_c6 = ophthalmic_id if ophthalmic_id is not None else c5_c6
        c6_c7 = p_com_a_id if p_com_a_id is not None else c6_c7

        # Adjust C5 if C6 has moved.
        c5_length_new = c5_c6 - c4_c5
        if c5_length_new != c5_length:
            c4_c5 = c5_c6 - c5_length

    # Add C1 if C2 is not at centerline ID=0
    if very_short_model:
        interfaces = {"C4": c3_c4, "C5": c4_c5, "C6": c5_c6, "C7": c6_c7}
    elif short_model:
        interfaces = {"C3": c2_c3, "C4": c3_c4, "C5": c4_c5, "C6": c5_c6, "C7": c6_c7}
    elif c1_c2 == 0:
        interfaces = {"C2": c1_c2, "C3": c2_c3, "C4": c3_c4, "C5": c4_c5, "C6": c5_c6, "C7": c6_c7}
    else:
        interfaces = {"C1": 0, "C2": c1_c2, "C3": c2_c3, "C4": c3_c4, "C5": c4_c5, "C6": c5_c6, "C7": c6_c7}

    # Set landmarks
    landmarks = {}
    for k, v in interfaces.items():
        landmarks[k] = line.GetPoint(int(v))

    # Map landmarks to initial centerline
    landmarks = map_landmarks(landmarks, centerline, algorithm="kjeldsberg")

    # Save landmarks
    print("-- Case was successfully landmarked.")
    print("-- Number of landmarks (Segments): %s" % len(landmarks))
    try:
        os.remove(base_path + "_landmark_kjeldsberg_vmtk.particles")
    except:
        pass

    if landmarks is not None:
        write_parameters(landmarks, base_path)
        create_particles(base_path, algorithm="kjeldsberg", method="vmtk")

    return landmarks


def adjust_c5_c6_interface(c4_c5, c5_c6, c6_c7, curvilinear_coordinates):
    """
    Adjust C5/C6 interface based on length requirement

    Args:
        c4_c5 (int): ID of C4/C5 interface
        c5_c6 (int): ID of C5/C6 interface
        c6_c7 (int): ID of C6/C7 interface
        curvilinear_coordinates (ndarray): Array of curvilinear coordinates

    Returns:
        c5_c6 (int): Adjusted C5/C6 interface
    """
    len_c6 = la.norm(curvilinear_coordinates[c6_c7] - curvilinear_coordinates[c5_c6])
    len_c5 = la.norm(curvilinear_coordinates[c5_c6] - curvilinear_coordinates[c4_c5])
    while len_c6 * 2 / 3 < len_c5:
        dx = c5_c6 - c4_c5
        c5_c6 -= int(0.1 * dx)
        len_c5 = la.norm(curvilinear_coordinates[c5_c6] - curvilinear_coordinates[c4_c5])
        len_c6 = la.norm(curvilinear_coordinates[c6_c7] - curvilinear_coordinates[c5_c6])
    return c5_c6


def adjust_c6_c7_interface(c5_c6, c6_c7, curvilinear_coordinates):
    """
    Adjust C6/C7 interface based on length requirement

    Args:
        c5_c6 (int): ID of C5/C6 interface
        c6_c7 (int): ID of C6/C7 interface
        curvilinear_coordinates (ndarray): Array of curvilinear coordinates

    Returns:
        c6_c7 (int): Adjusted C6/C7 interface
    """
    len_c6 = la.norm(curvilinear_coordinates[c6_c7] - curvilinear_coordinates[c5_c6])
    len_c7 = la.norm(curvilinear_coordinates[-1] - curvilinear_coordinates[c6_c7])
    while len_c6 < len_c7:
        dx = c6_c7 - c5_c6
        c6_c7 += int(0.1 * dx)
        len_c6 = la.norm(curvilinear_coordinates[c6_c7] - curvilinear_coordinates[c5_c6])
        len_c7 = la.norm(curvilinear_coordinates[-1] - curvilinear_coordinates[c6_c7])
    return c6_c7


def find_c6_c7_interface(c5_c6, max_point_ids, curvilinear_coordinates):
    """
    Find C6/C7 interface based on curvature maximum

    Args:
        c5_c6 (int): ID of C5/C6 interface
        max_point_ids (ndarray): Array with curvature maxima
        curvilinear_coordinates (ndarray): Array of curvilinear coordinates

    Returns:
        c6_c7 (int): C6/C7 interface
    """
    c6_c7_curvature_index = (max_point_ids > c5_c6).nonzero()[0]
    if c6_c7_curvature_index.size > 0:
        c6_c7 = max_point_ids[c6_c7_curvature_index.max()]
    else:
        c6_c7 = int((len(curvilinear_coordinates) + c5_c6) * 0.5)
    return c6_c7


def find_c5_c6_interface(c4_c5, tangent):
    """
    Find C5/C6 interface based on angle requirement (>30 degrees)

    Args:
        c4_c5 (int): ID of C4/C5 interface
        tangent (ndarray): Array of tangent vectors along centerline

    Returns:
        c5_c6 (int): C5/C6 interface
    """
    t45 = tangent[c4_c5]
    c5_c6 = np.asarray(int(c4_c5 * 1.02))
    t56 = tangent[c5_c6]
    alpha = np.arccos(abs(np.dot(t45, t56)) / (la.norm(t45) * la.norm(t56))) * 180 / np.pi
    alpha_tolerance = 30
    while alpha < alpha_tolerance:
        c5_c6 = np.array(int(c5_c6 * 1.02))
        t56 = tangent[c5_c6]
        alpha = np.arccos(abs(np.dot(t45, t56)) / (la.norm(t45) * la.norm(t56))) * 180 / np.pi
    return c5_c6


def adjust_c1_c2_interface(c1_c2, c2_c3, tangent):
    """
    Adjust C1/C2 interface based on angle requirement (>60 degrees)

    Args:
        c1_c2 (int): ID of C1/C2 interface
        c2_c3 (int): ID of C2/C3 interface
        tangent (ndarray): Array of tangent vectors along centerline

    Returns:
        c1_c2 (int): Adjusted C1/C2 interface
    """

    t12 = tangent[c1_c2]
    t23 = tangent[c2_c3]
    alpha = np.arccos(abs(np.dot(t12, t23)) / (la.norm(t12) * la.norm(t23))) * 180 / np.pi
    alpha_tolerance = 60
    c1_c2_min = 0
    while alpha < alpha_tolerance and c1_c2 > c1_c2_min:
        c1_c2 = np.array(int(c1_c2 * 0.95))
        t12 = tangent[c1_c2]
        alpha = np.arccos(abs(np.dot(t12, t23)) / (la.norm(t12) * la.norm(t23))) * 180 / np.pi
    if c1_c2 < 0:
        c1_c2 = np.asarray(0)
    return c1_c2


def find_c1_c2_interface(max_point_ids, min_point_ids, start):
    """
    Find C1/C2 interface based on curvature minimum

    Args:
        max_point_ids (ndarray): Array with curvature maxima
        min_point_ids (ndarray): Array with curvature mimuma
        start (int): ID to start searching from

    Returns:
        c1_c2 (int): C1/C2 interface
    """
    stop = start
    start_ids = (max_point_ids < stop).nonzero()[0]
    if start_ids.size > 0:
        start = max_point_ids[(max_point_ids < stop).nonzero()[0][-1]]
        c1_c2_index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]
        c1_c2 = min_point_ids[c1_c2_index[0]]
    else:
        c1_c2 = np.array(0)
    return c1_c2


def adjust_c2_c3_interface(c2_c3, c3_c4, tangent):
    """
    Adjust C2/C3 based on vertical interface criteria relative to C3/C4

    Args:
        c2_c3 (int): ID of C2/C3 interface
        c3_c4 (int): ID of C3/C4 interface
        tangent (ndarray): Array of tangent vectors along centerline

    Returns:
        c2_c3 (int): Adjusted C2/C3 interface
    """
    t23 = tangent[c2_c3]
    t34 = tangent[c3_c4]
    alpha = np.arccos(abs(np.dot(t23, t34)) / (la.norm(t23) * la.norm(t34))) * 180 / np.pi
    alpha_tolerance = 60
    while alpha < alpha_tolerance:
        c2_c3 = np.array(int(c2_c3 * 0.95))
        t23 = tangent[c2_c3]
        alpha = np.arccos(abs(np.dot(t23, t34)) / (la.norm(t23) * la.norm(t34))) * 180 / np.pi
    return c2_c3


def find_c2_c3_interface(max_point_ids, min_point_ids, start, short_model, very_short_model):
    """
    Find C2/C3 based on curvature minimum.

    Args:
        max_point_ids (ndarray): Array with curvature maxima
        min_point_ids (ndarray): Array with curvature mimuma
        start (int): ID to start searching from
        short_model (boolean): True if model does not contain C1 or C2
        very_short_model (boolean): True if model does not contain C1, C2 or C3

    Returns:
        c2_c3 (int): C2/C3 interface
        short_model (boolean): True if model does not contain C1 or C2
        start (int): ID to start searching from
    """
    stop = start

    if not very_short_model and (max_point_ids < stop).nonzero()[0].size > 0:
        start = max_point_ids[(max_point_ids < stop).nonzero()[0][-1]]
        c2_c3_index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]
        c2_c3 = min_point_ids[c2_c3_index[0]]
    else:
        short_model = True
        c2_c3 = 0
    return c2_c3, short_model, start


def find_c3_c4_interface(c4_c5, max_point_ids, min_point_ids, very_short_model):
    """
    Find C3/C4 based on curvature minimum between two local maxima.

    Args:
        c4_c5 (int): ID at C4/C5 interface
        max_point_ids (ndarray): Array with curvature maxima
        min_point_ids (ndarray): Array with curvature mimuma
        very_short_model (boolean): True if model does not contain C1, C2 or C3

    Returns:
        c3_c4 (int): C3/C4 interface
        start (int): ID to start searching from
        very_short_model (boolean): True if model does not contain C1, C2 or C3
    """
    anterior_bend_peak_id = max_point_ids[(max_point_ids <= c4_c5).nonzero()[0].max()]
    region_points = (max_point_ids < anterior_bend_peak_id).nonzero()[0]
    if region_points.size > 1:
        region_ids = region_points[-2:]
        start = max_point_ids[region_ids[0]]
        stop = max_point_ids[region_ids[1]]
    elif region_points.size == 1:
        start = np.asarray(0)
        stop = max_point_ids[region_points[0]]
    else:
        return 0, np.asarray(0), True

    c3_c4_index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]
    if c3_c4_index.size > 0:
        c3_c4 = min_point_ids[c3_c4_index[0]]
    else:
        c3_c4 = 0
        very_short_model = True

    return c3_c4, start, very_short_model


def find_c4_c5_interface(curvilinear_coordinates, coronal_coordinates):
    """
    Find C4/C5 at the maximum coronal coordinate along the centerline.

    Args:
        curvilinear_coordinates (ndarray): Array of curvilinear coordinates
        coronal_coordinates (str): Array of coronal coordinates
    Returns:
        c4_c5 (int): C4/C5 interface
    """
    value_index = coronal_coordinates[argrelextrema(coronal_coordinates, np.less_equal)[0]].min()
    max_coronal_bend_id = np.array(coronal_coordinates.tolist().index(value_index))
    search_tolerance = curvilinear_coordinates[-1] * 3 / 4
    if abs(curvilinear_coordinates[max_coronal_bend_id] - curvilinear_coordinates[-1]) > search_tolerance:
        print("-- Sanity check failed, checking for maximums")

        value_index = coronal_coordinates[argrelextrema(coronal_coordinates, np.greater_equal)[0]].max()
        max_coronal_bend_id = np.array(coronal_coordinates.tolist().index(value_index))
        if abs(curvilinear_coordinates[max_coronal_bend_id] - curvilinear_coordinates[-1]) > search_tolerance:
            print("-- Sanity check failed, no anterior bend in model")
            sys.exit(1)

    c4_c5 = max_coronal_bend_id
    return c4_c5

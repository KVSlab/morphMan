##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

## This software is distributed WITHOUT ANY WARRANTY; without even
## the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE. See the above copyright notices for more information.

import os

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.signal import argrelextrema

from morphman.automated_landmarking.automated_landmarking_tools import (
    create_particles,
    map_landmarks,
    spline_centerline_and_compute_geometric_features,
)
from morphman.common.centerline_operations import get_curvilinear_coordinate
from morphman.common.tools_common import write_parameters
from morphman.common.vmtk_wrapper import (
    vmtk_compute_geometric_features,
    vmtk_resample_centerline,
)
from morphman.common.vtk_wrapper import get_point_data_array


def landmarking_piccinelli(
    centerline,
    base_path,
    approximation_method,
    algorithm,
    resampling_step,
    smooth_line,
    nknots,
    smoothing_factor_curv,
    smoothing_factor_torsion,
    iterations,
):
    """
    Perform landmarking of an input centerline to
    identify different segments along the vessel.
    Landmarking algorithm based on Bogunovic et al. (2012).
    Uses curvature and torsion to objectively subdivide
    the siphon into bends.
    Subdivision of individual siphon bends is
    performed by identifying locations of curvature and torsion peaks along
    the siphon and defining a bend for each curvature peak delimited by 2
    enclosing (proximal and distal) torsion peaks.

    Args:
        centerline (vtkPolyData): Centerline data points.
        base_path (str): Location of case to landmark.
        approximation_method (str): Method used for computing curvature.
        algorithm (str): Name of landmarking algorithm.
        resampling_step (float): Resampling step. Is None if no resampling.
        smooth_line (bool): Smooths centerline with VMTK if True.
        nknots (int): Number of knots for B-splines.
        smoothing_factor_curv (float): Smoothing factor for computing curvature.
        smoothing_factor_torsion (float): Smoothing factor for computing torsion.
        iterations (int): Number of smoothing iterations.

    Returns:
        landmarks (dict): Landmarking interfaces as points.
    """

    if resampling_step is not None:
        centerline = vmtk_resample_centerline(centerline, resampling_step)

    if approximation_method == "spline":
        line, max_point_curv_ids, _ = spline_centerline_and_compute_geometric_features(
            centerline, smooth_line, nknots
        )

        # Get curvature and torsion, find peaks
        curvature = get_point_data_array("Curvature", line)
        torsion = get_point_data_array("Torsion", line)
        torsion_smooth = gaussian_filter(torsion, 10)
        max_point_tor_ids = list(argrelextrema(abs(torsion_smooth), np.greater)[0])

    elif approximation_method == "vmtk":
        line = centerline
        line_curv = vmtk_compute_geometric_features(
            centerline,
            True,
            output_smoothed=False,
            factor=smoothing_factor_curv,
            iterations=iterations,
        )
        line_tor = vmtk_compute_geometric_features(
            centerline,
            True,
            output_smoothed=False,
            factor=smoothing_factor_torsion,
            iterations=iterations,
        )
        # Get curvature and torsion, find peaks
        curvature = get_point_data_array("Curvature", line_curv)
        torsion = get_point_data_array("Torsion", line_tor)

        # Smooth torsion curve to remove noise
        torsion_smooth = gaussian_filter(torsion, 25)

        # Find maximum curvature and torsion
        max_point_curv_ids = list(argrelextrema(curvature, np.greater)[0])
        max_point_tor_ids = list(argrelextrema(abs(torsion_smooth), np.greater)[0])

    else:
        raise ValueError(
            "ERROR: Selected method for computing curvature / torsion not available"
            + "\nPlease select between 'spline' and 'vmtk'"
        )

    # Extract local curvature minimums
    length = get_curvilinear_coordinate(line)

    # Ignore points close to the ends of the siphon
    for i in max_point_curv_ids:
        if length[i] in length[-10:] or length[i] in length[:10]:
            max_point_curv_ids.remove(i)

    # Remove curvature and torsion peaks too close to each other
    tolerance = 70
    dist = []
    dist_tor = []
    for i in range(len(max_point_curv_ids) - 1):
        dist.append(max_point_curv_ids[i + 1] - max_point_curv_ids[i])
    for i in range(len(max_point_tor_ids) - 1):
        dist_tor.append(max_point_tor_ids[i + 1] - max_point_tor_ids[i])

    # Remove curvature maxima which are saddle points within a bend
    curv_remove_ids = []
    for i, dx in enumerate(dist):
        if dx < tolerance:
            curv1 = curvature[max_point_curv_ids[i]]
            curv2 = curvature[max_point_curv_ids[i + 1]]
            if curv1 > curv2:
                curv_remove_ids.append(max_point_curv_ids[i + 1])
            else:
                curv_remove_ids.append(max_point_curv_ids[i])

    # Remove torsion maxima which are saddle points within a bend
    tor_remove_ids = []
    for i, dx in enumerate(dist_tor):
        if dx < tolerance:
            tor1 = torsion_smooth[max_point_tor_ids[i]]
            tor2 = torsion_smooth[max_point_tor_ids[i + 1]]
            if tor1 > tor2:
                tor_remove_ids.append(max_point_tor_ids[i + 1])
            else:
                tor_remove_ids.append(max_point_tor_ids[i])

    # Filter out curvature and torsion maxima saddle points
    max_point_curv_ids = [ID for ID in max_point_curv_ids if ID not in curv_remove_ids]
    max_point_tor_ids = [ID for ID in max_point_tor_ids if ID not in tor_remove_ids]

    # Compute and extract interface points
    interfaces = find_interface(max_point_curv_ids, max_point_tor_ids)
    landmarks = {}
    for k, v in interfaces.items():
        landmarks[k] = line.GetPoints().GetPoint(int(v))

    # Map landmarks to initial centerline
    landmarks = map_landmarks(landmarks, centerline, algorithm)

    # Save landmarks
    print("-- Case was successfully landmarked.")
    print("-- Number of landmarks (Segments): %s" % len(landmarks))
    try:
        os.remove(
            base_path + "_landmark_piccinelli_%s.particles" % approximation_method
        )
    except FileNotFoundError:
        pass
    if landmarks is not None:
        write_parameters(landmarks, base_path)
        create_particles(base_path, algorithm, approximation_method)

    return landmarks


def find_interface(max_point_curv_ids, max_point_tor_ids):
    """
    Find interfaces between bends as defined by Piccinelli et al.
    A bend is defined by a curvature maximum
    bounded by two enclosing torsion maxima.

    Args:
        max_point_curv_ids (ndarray): Array of curvature maximum
        max_point_tor_ids (ndarray): Array of torsion maximum

    Returns:
        interface (dict): Dictionary of interfaces between bends
    """
    found = False
    interface = {}
    k = 0
    start_id = 0
    for c in max_point_curv_ids:
        for i in range(start_id, len(max_point_tor_ids) - 1):
            if max_point_tor_ids[i] < c < max_point_tor_ids[i + 1] and not found:
                interface["bend%s" % (k + 1)] = np.array([max_point_tor_ids[i]])
                k += 1
                interface["bend%s" % (k + 1)] = np.array([max_point_tor_ids[i + 1]])
                k += 1
                start_id = i + 1
                found = True
        found = False

    return interface

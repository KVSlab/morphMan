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
    get_centerline_coordinates,
    get_maximum_coronal_coordinate,
    map_landmarks,
    spline_centerline_and_compute_geometric_features,
)
from morphman.common.centerline_operations import (
    compute_discrete_derivatives,
    get_curvilinear_coordinate,
    get_k1k2_basis,
)
from morphman.common.tools_common import write_parameters
from morphman.common.vmtk_wrapper import (
    vmtk_compute_centerline_attributes,
    vmtk_compute_geometric_features,
    vmtk_resample_centerline,
)
from morphman.common.vtk_wrapper import get_point_data_array


def landmarking_bogunovic(
    centerline,
    base_path,
    approximation_method,
    algorithm,
    resampling_step,
    smooth_line,
    nknots,
    smoothing_factor,
    iterations,
    coronal_axis,
):
    """
    Perform landmarking of an input centerline to
    identify different segments along the vessel.
    Landmarking algorithm based on Bogunovic (2012).
    Bends are identified where the angle in the k1-k2 basis
    exceeds a certain limit.

    Args:
        centerline (vtkPolyData): Centerline data points.
        base_path (str): Location of case to landmark.
        approximation_method (str): Method used for computing curvature.
        algorithm (str): Name of landmarking algorithm.
        resampling_step (float): Resampling step. Is None if no resampling.
        smooth_line (bool): Smooths centerline with VMTK if True.
        nknots (int): Number of knots for B-splines.
        smoothing_factor (float): Smoothing factor used in VMTK
        iterations (int): Number of smoothing iterations.
        coronal_axis (str) : Axis determining coronal coordinate

    Returns:
        landmarks (dict): Landmarking interfaces as points.
    """

    if resampling_step is not None:
        centerline = vmtk_resample_centerline(centerline, length=resampling_step)

    if approximation_method == "vmtk":
        line = vmtk_compute_centerline_attributes(centerline)
        line = vmtk_compute_geometric_features(
            line, smooth_line, factor=smoothing_factor, iterations=iterations
        )
        curvature = get_point_data_array("Curvature", line)
        curvature = gaussian_filter(curvature, 2)

    elif approximation_method == "disc":
        neigh = 20
        line = vmtk_compute_centerline_attributes(centerline)
        line = vmtk_compute_geometric_features(
            line, smooth_line, factor=smoothing_factor, iterations=iterations
        )
        line, curvature__ = compute_discrete_derivatives(line, neigh=neigh)
        curvature = []
        for c in curvature__:
            curvature.append([c])
        curvature = np.array(curvature)

    elif approximation_method == "spline":
        line, max_point_ids, min_point_ids = (
            spline_centerline_and_compute_geometric_features(
                centerline, smooth_line, nknots
            )
        )
        curvature = get_point_data_array("Curvature", line)

    if approximation_method != "spline":
        max_point_ids = list(argrelextrema(curvature, np.greater)[0])
        min_point_ids = list(argrelextrema(curvature, np.less)[0])
        get_k1k2_basis(curvature, line)

    curvilinear_coordinates = get_curvilinear_coordinate(line)
    k1 = get_point_data_array("k1", line)
    k2 = get_point_data_array("k2", line)

    # Remove a min / max point that is in reality a saddle point
    for i in min_point_ids:
        for j in max_point_ids:
            if abs(i - j) < 5 and abs(curvature[i] - curvature[j]) < 0.01:
                if i in min_point_ids:
                    min_point_ids.remove(i)
                if j in max_point_ids:
                    max_point_ids.remove(j)

    # Find local curvature max in the k1-k2 space, and compute angle between the points
    k1_points = k1[max_point_ids]
    k2_points = k2[max_point_ids]
    k_points = np.zeros((k1_points.shape[0], 2))
    k_points[:, 0] = k1_points[:, 0]
    k_points[:, 1] = k2_points[:, 0]
    theta = np.zeros(k1_points.shape[0] - 1)
    for i in range(theta.shape[0]):
        a = k_points[i, :] / np.sqrt(np.sum(k_points[i, :] ** 2))
        b = k_points[i + 1, :] / np.sqrt(np.sum(k_points[i + 1, :] ** 2))
        theta[i] = np.arccos(np.dot(a, b))
        theta[i] = theta[i] * 180 / np.pi

    # Tolerance parameters from Bogunovic et al. (2012)
    tol_anterior_posterior = 60
    tol_superior_anterior = 45
    tol_posterior_inferior = 45
    tol_inferior_end = 110

    # Find coronal coordinate and maximum (within anterior bend)
    coordinates = get_centerline_coordinates(line, curvilinear_coordinates)
    coronal_coordinate = coordinates[coronal_axis]
    max_coronal_coordinate_id = get_maximum_coronal_coordinate(
        coronal_coordinate, curvilinear_coordinates
    )

    interfaces = {}
    min_point_ids = np.array(min_point_ids)

    # Find interfaces between bends based on angles between curvature vectors and tolerance parameters
    index = np.array((max_coronal_coordinate_id > max_point_ids).nonzero()[0]).max()
    start = find_interface(
        index,
        -1,
        tol_anterior_posterior,
        "anterior_posterior",
        theta,
        max_point_ids,
        min_point_ids,
        interfaces,
    )
    if start is None:
        return None
    start = find_interface(
        start,
        -1,
        tol_posterior_inferior,
        "posterior_inferior",
        theta,
        max_point_ids,
        min_point_ids,
        interfaces,
    )
    if start is None:
        return None
    find_interface(
        start,
        -1,
        tol_inferior_end,
        "inferior_end",
        theta,
        max_point_ids,
        min_point_ids,
        interfaces,
    )
    start = find_interface(
        index + 1,
        1,
        tol_superior_anterior,
        "superior_anterior",
        theta,
        max_point_ids,
        min_point_ids,
        interfaces,
    )
    if start is None:
        return None

    # Set landmarks following Bogunovic
    landmarks = {}
    for k, v in interfaces.items():
        landmarks[k] = line.GetPoint(int(v))

    # Map landmarks to initial centerline
    landmarks = map_landmarks(landmarks, centerline, algorithm)

    # Save landmarks
    print("-- Case was successfully landmarked.")
    print("-- Number of landmarks (Segments): %s" % len(landmarks))
    try:
        os.remove(base_path + "_landmark_bogunovic_%s.particles" % approximation_method)
    except FileNotFoundError:
        pass

    if landmarks is not None:
        write_parameters(landmarks, base_path)
        create_particles(base_path, algorithm, approximation_method)

    return landmarks


def find_interface(
    start, direction, tol, part, theta, max_point_ids, min_point_ids, interfaces
):
    """
    Find landmark interfaces according to Bogunovic.
    Based on local curvature extrema computed in the k1-k2 space.
    Angles between curvature vectors determine interfaces based on a pre-defined tolerance.

    Args:
        start (int): Index where search for curvature max starts
        direction (int): Direction to search (1 or -1)
        tol (double): Tolerance on angle between curvature vectors
        part (str): Name of interface
        theta (ndarray): Array of angles between curvature vectors
        max_point_ids (ndarray): Local curvature maximum
        min_point_ids (ndarray): Local curvature minimum
        interfaces (dict): Dictionary of landmarks/interfaces and coresponding points

    Returns:
        i (int): Index of curvature maximum in bend of interest
    """
    stop = direction if direction == -1 else theta.shape[0]
    success = False
    for i in range(start - 1, stop, direction):
        if theta[i] > tol:
            success = True
            break

    if success:
        start = max_point_ids[i]
        stop = max_point_ids[i + 1]
        index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0].max()
        min_point = min_point_ids[index]
        interfaces[part] = min_point

    elif not success and part == "superior_anterior":
        print(
            "-- Where not able to identify the interface between the "
            + "anterior and superior bend. Check the coronal coordinates"
        )
        return None

    elif not success and part != "inferior_end":
        print(
            "-- The geometry is to short to be classified with superior"
            + ", anterior, posterior and inferior."
        )
        return None

    elif not success and part == "inferior_end":
        interfaces["inf_end"] = np.array([0])
        i = 0
        print(
            "-- End of inferior is at the end of the geometry, this might"
            + "affect the geometry stats"
        )
    else:
        print("-- Something happened, idea: some bend ended at the last point")
        return None

    return i

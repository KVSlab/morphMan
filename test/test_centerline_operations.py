##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

# noinspection PyUnresolvedReferences
from IPython import embed

from .fixtures import surface_paths
from morphman.common.centerline_operations import *


def test_diverging_point_on_centerline(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline1 = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    centerline2 = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 1)

    len1 = centerline1.GetNumberOfPoints()
    len2 = centerline2.GetNumberOfPoints()

    tol = 1E-5

    # Sanity check
    index = get_diverging_point_id(centerline1, centerline1, tol)
    assert index == (len1 - 1)

    index = get_diverging_point_id(centerline1, centerline2, tol)
    assert index < len1 and index < len2


def test_curvilinear_coordinates(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    coordinates = get_curvilinear_coordinate(centerline)

    length = len(coordinates)

    assert length == centerline.GetNumberOfPoints()


def test_centerline_tolerance(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    tolerance = get_centerline_tolerance(centerline)

    diff = np.asarray(centerline.GetPoint(10)) - np.asarray(centerline.GetPoint(9))
    tolerance_max = np.linalg.norm(diff)

    assert tolerance < tolerance_max


def test_clip_diverging_centerline(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    n_points = centerline.GetNumberOfPoints()
    start = int(0.9 * n_points)
    p0 = centerline.GetPoint(start)

    clipped_centerline = get_clipped_diverging_centerline(centerline, p0, n_points - 1)

    n_clipped_points = clipped_centerline.GetNumberOfPoints()

    diff = abs(n_points * 0.9 - n_clipped_points)

    assert diff < 5


def test_get_first_line_to_change(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")
    full_line_to_change = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    line_to_change, _, _, _, _ = get_line_to_change(None, centerlines, "first_line", None, None, None)

    n_points_full = full_line_to_change.GetNumberOfPoints()
    n_points = line_to_change.GetNumberOfPoints()

    p0 = np.asarray(line_to_change.GetPoint(0))
    p1 = np.asarray(line_to_change.GetPoint(0))

    pA = np.asarray(full_line_to_change.GetPoint(0))
    pB = np.asarray(full_line_to_change.GetPoint(0))

    tol = 1E-1

    assert n_points < n_points_full

    assert np.linalg.norm(p0 - pA) < tol
    assert np.linalg.norm(p1 - pB) < tol


def test_get_region_of_interest(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")
    n_points = centerlines.GetNumberOfPoints()

    start_id = int(n_points * 0.1)
    stop_id = int(n_points * 0.2)

    start = np.asarray(centerlines.GetPoint(start_id))
    stop = np.asarray(centerlines.GetPoint(stop_id))

    region_points = [start, stop]

    remaining_centerlines, diverging_centerlines, _, _, diverging_ids \
        = get_region_of_interest_and_diverging_centerlines(centerlines, region_points)

    assert centerlines.GetNumberOfLines() == 7
    assert remaining_centerlines.GetNumberOfLines() == 5
    assert diverging_centerlines.GetNumberOfLines() == 2

    assert len(diverging_ids) == 2


def test_discrete_derivatives(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    _, curvature = compute_discrete_derivatives(centerline)

    # Check if curvature computed for each point on centerline
    assert len(curvature) == centerline.GetNumberOfPoints()

    # Check no diverging values
    assert np.mean(curvature) < 0.5


def test_spline_centerline(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    _, curvature = compute_splined_centerline(
        centerline, get_curv=True, isline=True, nknots=15, get_stats=False, get_misr=False)

    # Check if curvature computed for each point on centerline
    assert len(curvature) == centerline.GetNumberOfPoints()

    # Check no diverging values
    assert np.mean(curvature) < 0.5


def test_get_end_point(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    last_point_true = np.asarray(centerline.GetPoint(centerline.GetNumberOfPoints() - 1))

    last_point = np.asarray(get_end_point(centerline))

    diff = np.linalg.norm(last_point - last_point_true)
    tol = 1E-16

    assert diff < tol

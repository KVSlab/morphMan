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



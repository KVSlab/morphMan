##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import numpy as np
import pytest

from .fixtures import surface_paths
from morphman import manipulate_area, read_command_line_area
from morphman.common.centerline_operations import get_curvilinear_coordinate
from morphman.common.common import get_path_names
from morphman.common.surface_operations import read_polydata, vmtk_compute_centerline_sections, get_point_data_array, \
    extract_single_line

"""
def test_area_linear(surface_paths):
    # Get default input
    common_input = read_command_line_area(surface_paths[0], surface_paths[1])

    # Set region points
    base_path = get_path_names(common_input['input_filepath'])
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    n = centerline.GetNumberOfPoints()
    region_points = list(centerline.GetPoint(int(n * 0.3))) + list(centerline.GetPoint(int(n * 0.4)))

    # Change default input
    common_input.update(dict(region_of_interest="commandline",
                             region_points=region_points,
                             method="linear",
                             smooth=False,
                             size=2,
                             percentage=50))

    # Run area variation
    manipulate_area(**common_input)

    # Set file paths
    base_path = get_path_names(common_input['input_filepath'])
    centerline_spline_path = base_path + "_centerline_spline.vtp"
    new_surface_path = common_input["output_filepath"]

    # Read data, and get new area
    surface = read_polydata(new_surface_path)
    centerline_spline = read_polydata(centerline_spline_path)
    new_centerline_area, _ = vmtk_compute_centerline_sections(surface,
                                                              centerline_spline)

    length = get_curvilinear_coordinate(centerline_spline)
    new_area = get_point_data_array("CenterlineSectionArea", new_centerline_area)
    linear_change = new_area[0] + (new_area[-1] - new_area[0]) * (length / length.max())

    # Check if the new area is within 2 % of the expected value
    assert (np.abs(new_area[:, 0] - linear_change) / linear_change).max() < 0.02


@pytest.mark.parametrize("ratio", [1.5, 3.0])
def test_area_variation(ratio, surface_paths):
    # Get default input
    common_input = read_command_line_area(surface_paths[0], surface_paths[1])

    # Set region points
    base_path = get_path_names(common_input['input_filepath'])
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    n = centerline.GetNumberOfPoints()
    region_points = list(centerline.GetPoint(int(n * 0.05))) + list(centerline.GetPoint(int(n * 0.5)))

    # Change default input
    common_input.update(dict(method="variation",
                             smooth=False,
                             region_of_interest="commandline",
                             region_points=region_points,
                             ratio=ratio,
                             beta=None))

    # Run area variation
    manipulate_area(**common_input)

    # Set file paths
    base_path = get_path_names(common_input['input_filepath'])
    centerline_spline_path = base_path + "_centerline_spline.vtp"
    new_surface_path = common_input["output_filepath"]

    # Read data, and get new area
    surface = read_polydata(new_surface_path)
    centerline_spline = read_polydata(centerline_spline_path)
    new_centerline_area, _ = vmtk_compute_centerline_sections(surface,
                                                              centerline_spline)

    new_area = get_point_data_array("CenterlineSectionArea", new_centerline_area)

    # Check if the new ratio holds
    assert abs(ratio - new_area.max() / new_area.min()) < 0.15


def test_create_stenosis(surface_paths):
    # Get default input
    common_input = read_command_line_area(surface_paths[0], surface_paths[1])

    # Get region points
    base_path = get_path_names(common_input['input_filepath'])
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    n = centerline.GetNumberOfPoints()
    region_point = list(centerline.GetPoint(int(n * 0.4)))

    # Change default input
    common_input.update(dict(region_of_interest="commandline",
                             region_points=region_point,
                             method="stenosis",
                             size=1.0,
                             percentage=50))

    # Create a stenosis
    manipulate_area(**common_input)

    # Import old area and splined centerline for region of interest
    base_path = get_path_names(common_input['input_filepath'])
    centerline_spline_path = base_path + "_centerline_spline.vtp"
    centerline_area_spline_path = base_path + "_centerline_area_spline.vtp"
    new_surface_path = common_input["output_filepath"]
    surface = read_polydata(new_surface_path)
    centerline_area = read_polydata(centerline_area_spline_path)
    centerline_spline = read_polydata(centerline_spline_path)
    new_centerline_area, _ = vmtk_compute_centerline_sections(surface,
                                                              centerline_spline)
    old_area = get_point_data_array("CenterlineSectionArea", centerline_area)
    new_area = get_point_data_array("CenterlineSectionArea", new_centerline_area)

    # Check if there is a 50 % narrowing
    assert abs((np.sqrt(new_area / old_area)).min() - 0.5) < 0.05


"""


def test_create_bulge(surface_paths):
    # Get default input
    common_input = read_command_line_area(surface_paths[0], surface_paths[1])

    # Get region points
    base_path = get_path_names(common_input['input_filepath'])
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    n = centerline.GetNumberOfPoints()
    region_point = list(centerline.GetPoint(int(n * 0.4)))

    # Set problem specific parameters
    common_input.update(dict(method="bulge",
                             region_of_interest="commandline",
                             region_points=region_point,
                             size=1.0,
                             percentage=50,
                             smooth=False))

    # Create a stenosis
    manipulate_area(**common_input)

    # Import old area and splined centerline for region of interest
    base_path = get_path_names(common_input['input_filepath'])
    centerline_spline_path = base_path + "_centerline_spline.vtp"
    centerline_area_spline_path = base_path + "_centerline_area_spline.vtp"
    new_surface_path = common_input["output_filepath"]
    surface = read_polydata(new_surface_path)
    centerline_area = read_polydata(centerline_area_spline_path)
    centerline_spline = read_polydata(centerline_spline_path)
    new_centerline_area, _ = vmtk_compute_centerline_sections(surface,
                                                              centerline_spline)
    old_area = get_point_data_array("CenterlineSectionArea", centerline_area)
    new_area = get_point_data_array("CenterlineSectionArea", new_centerline_area)

    # Check if there is a 50 % narrowing
    assert abs((np.sqrt(new_area / old_area)).max() - 1.5) < 0.05


@pytest.mark.parametrize("percentage", [-15, 15])
def test_inflation_and_deflation_of_area(surface_paths, percentage):
    # Get default input
    common_input = read_command_line_area(surface_paths[0], surface_paths[1])

    # Change the default input
    common_input.update(dict(method="area",
                             region_of_interest="first_line",
                             percentage=percentage,
                             smooth=False))

    # Perform area manipulation
    manipulate_area(**common_input)

    # Import old area and splined centerline for region of interest
    base_path = get_path_names(common_input['input_filepath'])
    centerline_spline_path = base_path + "_centerline_spline.vtp"
    centerline_area_spline_path = base_path + "_centerline_area_spline.vtp"
    surface = read_polydata(common_input["output_filepath"])
    centerline_area = read_polydata(centerline_area_spline_path)
    centerline_spline = read_polydata(centerline_spline_path)
    new_centerline_area, _ = vmtk_compute_centerline_sections(surface,
                                                              centerline_spline)
    old_area = get_point_data_array("CenterlineSectionArea", centerline_area)
    new_area = get_point_data_array("CenterlineSectionArea", new_centerline_area)

    # Exclude first 5 %
    old_area = old_area[int(0.1 * old_area.shape[0]):]
    new_area = new_area[int(0.1 * new_area.shape[0]):]

    # Change in radius
    ratio = np.sqrt(new_area / old_area)

    # Check if the altered area is equal has change according to percentage
    assert np.mean(np.abs(ratio - (1 + percentage * 0.01))) < 0.05

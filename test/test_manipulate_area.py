##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import pytest
import numpy as np

from .fixtures import common_input
from morphman.common import read_polydata, vmtk_compute_centerline_sections, get_array, \
                            get_path_names, extract_single_line, get_curvilinear_coordinate
from morphman import manipulate_area

def test_area_linear(common_input):
    base_path = get_path_names(common_input['input_filepath'])
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    n = centerline.GetNumberOfPoints()
    region_points = list(centerline.GetPoint(int(n * 0.3))) + list(centerline.GetPoint(int(n * 0.4)))
    common_input['region_of_interest'] = "commandline"
    common_input['region_points'] = region_points
    common_input.update(dict(method = "linear",
                             stenosis_length = 2,
                             percentage = 50,
                             ratio = None,
                             beta = None))

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
    new_area = get_array("CenterlineSectionArea", new_centerline_area)
    linear_change = new_area[0] + (new_area[-1] - new_area[0]) * (length / length.max())

    #linear_change.dump("linear_change.np")
    #new_area.dump("area.np")

    # Check if the new ratio holds
    assert np.abs(new_area[:,0] - linear_change).max() < 0.1


@pytest.mark.parametrize("ratio", [1.5, 3.0])
def test_area_variation(ratio, common_input):
    common_input.update(dict(method = "variation",
                             region_points = None,               # Inactive
                             region_of_interest = "first_line",
                             stenosis_length = 0,                # Inactive
                             percentage = 0,                     # Inactive
                             ratio = ratio,
                             beta = None))

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

    new_area = get_array("CenterlineSectionArea", new_centerline_area)

    # Check if the new ratio holds
    assert abs(ratio - new_area.max() / new_area.min()) < 0.1


def test_create_stenosis(common_input):
    # Get region points
    base_path = get_path_names(common_input['input_filepath'])
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    n = centerline.GetNumberOfPoints()
    region_point = list(centerline.GetPoint(int(n * 0.4)))
    common_input['region_of_interest'] = "commandline"
    common_input['region_points'] = region_point

    # Set problem specific parameters
    common_input.update(dict(method = "stenosis",
                             stenosis_length = 1.0,
                             percentage = 50,
                             ratio = None,   # Inactive
                             beta = None))   # Inactive

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
    old_area = get_array("CenterlineSectionArea", centerline_area)
    new_area = get_array("CenterlineSectionArea", new_centerline_area)

    # Check if there is a 50 % narrowing
    assert abs((np.sqrt(new_area / old_area)).min() - 0.5) < 0.05


def test_create_bulge(common_input):
    # Get region points
    base_path = get_path_names(common_input['input_filepath'])
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    n = centerline.GetNumberOfPoints()
    region_point = list(centerline.GetPoint(int(n * 0.4)))
    common_input['region_of_interest'] = "commandline"
    common_input['region_points'] = region_point

    # Set problem specific parameters
    common_input.update(dict(method = "bulge",
                             stenosis_length = 1.0,
                             percentage = 50,
                             ratio = None,   # Inactive
                             beta = None))   # Inactive

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
    old_area = get_array("CenterlineSectionArea", centerline_area)
    new_area = get_array("CenterlineSectionArea", new_centerline_area)

    # Check if there is a 50 % narrowing
    assert abs((np.sqrt(new_area / old_area)).max() - 1.5) < 0.05



@pytest.mark.parametrize("percentage",[-15, 15])
def test_inflation_and_deflation_of_area(common_input, percentage):
    common_input.update(dict(method = "area",
                             region_points = None,               # Inactive
                             region_of_interest = "first_line",
                             stenosis_length = 0,                # Inactive
                             percentage = percentage,
                             ratio = None,                       # Inactive
                             beta = None))                       # Inactive

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
    old_area = get_array("CenterlineSectionArea", centerline_area)
    new_area = get_array("CenterlineSectionArea", new_centerline_area)

    # Exclude first 5 %
    old_area = old_area[int(0.1*old_area.shape[0]):]
    new_area = new_area[int(0.1*new_area.shape[0]):]

    # Change in radius
    ratio = np.sqrt(new_area / old_area)

    # Check if the altered area is equal has change according to percentage
    assert np.mean(np.abs(ratio - (1 + percentage * 0.01))) < 0.05

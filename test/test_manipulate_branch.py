##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

import numpy as np

from .fixtures import surface_paths
from morphman import manipulate_branch, read_command_line_branch, get_end_point
from morphman.common.centerline_operations import extract_single_line
from morphman.common.common import get_path_names
from morphman.common.surface_operations import read_polydata


def test_manipulate_branch_translation(surface_paths):
    # Get default input
    common_input = read_command_line_branch(surface_paths[0], surface_paths[1])

    # Opthalmic artery in model C0001
    branch_number = 2

    # No rotation around new surface normal vector
    azimuth_angle = 0

    # No rotation around new surface 'tangent' vector
    polar_angle = 0

    # New location of branch
    new_branch_location = (45.7, 37.6, 42.5)

    # Branch translation method
    translation_method = 'commandline'

    # Change default input
    common_input.update(
        dict(polar_angle=polar_angle, branch_to_manipulate_number=branch_number, branch_location=new_branch_location,
             translation_method=translation_method, azimuth_angle=azimuth_angle))

    # Run manipulate branch
    manipulate_branch(**common_input)

    # Set file paths
    base_path = get_path_names(common_input['input_filepath'])
    old_centerlines_path = base_path + "_centerline.vtp"
    new_centerlines_path = base_path + "_centerline_moved.vtp"

    # Read data, and get new area
    old_centerlines = read_polydata(old_centerlines_path)
    new_centerlines = read_polydata(new_centerlines_path)
    old_centerlines = [extract_single_line(old_centerlines, i) for i in range(old_centerlines.GetNumberOfLines())]
    new_centerlines = [extract_single_line(new_centerlines, i) for i in range(new_centerlines.GetNumberOfLines())]

    untouched_centerlines = compare_end_points(old_centerlines, new_centerlines)

    assert len(untouched_centerlines) < len(new_centerlines)


def test_manipulate_branch_translation_and_polar_rotation(surface_paths):
    # Get default input
    common_input = read_command_line_branch(surface_paths[0], surface_paths[1])

    # Arbitrary branch in model C0001
    branch_number = 0

    # No rotation around new surface normal vector
    azimuth_angle0 = 0

    # New location of branch
    new_branch_location = (47.0, 27.8, 54.4)

    # Branch translation method
    translation_method = 'commandline'

    # Change default input
    common_input.update(
        dict(branch_to_manipulate_number=branch_number, branch_location=new_branch_location,
             translation_method=translation_method))

    # Run without rotation around surface normal
    common_input.update(dict(azimuth_angle=azimuth_angle0))

    # Run area variation
    manipulate_branch(**common_input)

    # Set file paths
    base_path = get_path_names(common_input['input_filepath'])
    old_centerlines_path = base_path + "_centerline.vtp"
    new_centerlines_path = base_path + "_centerline_moved.vtp"

    # Read centerlines
    old_centerlines = read_polydata(old_centerlines_path)
    old_centerlines = [extract_single_line(old_centerlines, i) for i in range(old_centerlines.GetNumberOfLines())]

    new_centerlines0 = read_polydata(new_centerlines_path)
    new_centerlines0 = [extract_single_line(new_centerlines0, i) for i in range(new_centerlines0.GetNumberOfLines())]

    # Perform same manipulation WITH rotation around new surface normal.
    azimuth_angle1 = np.pi
    common_input.update(dict(azimuth_angle=azimuth_angle1))

    # Run manipulate branch
    manipulate_branch(**common_input)

    # Read centerlines
    new_centerlines_path = base_path + "_centerline_moved_and_rotated.vtp"
    new_centerlines1 = read_polydata(new_centerlines_path)
    new_centerlines1 = [extract_single_line(new_centerlines1, i) for i in range(new_centerlines1.GetNumberOfLines())]

    # Compare end point of all centerlines
    unchanged_centerlines = compare_end_points(new_centerlines0, new_centerlines1)
    unchanged_centerlines0 = compare_end_points(new_centerlines0, old_centerlines)
    unchanged_centerlines1 = compare_end_points(new_centerlines1, old_centerlines)

    assert len(unchanged_centerlines) < len(old_centerlines)
    assert len(unchanged_centerlines0) < len(old_centerlines)
    assert len(unchanged_centerlines1) < len(old_centerlines)


def test_manipulate_branch_polar_rotation(surface_paths):
    # Get default input
    common_input = read_command_line_branch(surface_paths[0], surface_paths[1])

    # Arbitrary branch in model C0001
    branch_number = 2

    # No rotation around new surface normal vector
    polar_angle = np.pi / 4

    # Branch translation method
    translation_method = 'no_translation'

    # Change default input
    common_input.update(
        dict(branch_to_manipulate_number=branch_number, translation_method=translation_method))

    # Run with only rotation around original surface normal
    common_input.update(dict(polar_angle=polar_angle))

    # Run manipulate branch
    manipulate_branch(**common_input)

    # Set file paths
    base_path = get_path_names(common_input['input_filepath'])
    old_centerlines_path = base_path + "_centerline.vtp"
    new_centerlines_path = base_path + "_centerline_rotated.vtp"

    # Read centerlines
    old_centerlines = read_polydata(old_centerlines_path)
    old_centerlines = [extract_single_line(old_centerlines, i) for i in range(old_centerlines.GetNumberOfLines())]

    new_centerlines = read_polydata(new_centerlines_path)
    new_centerlines = [extract_single_line(new_centerlines, i) for i in range(new_centerlines.GetNumberOfLines())]

    # Compare end point of all centerlines
    unchanged_centerlines = compare_end_points(new_centerlines, old_centerlines)

    assert len(unchanged_centerlines) < len(old_centerlines)


def test_manipulate_branch_azimuthal_rotation(surface_paths):
    # Get default input
    common_input = read_command_line_branch(surface_paths[0], surface_paths[1])

    # Arbitrary branch in model C0001
    branch_number = 2

    # No rotation around new surface normal vector
    azimuth_angle = np.pi

    # Branch translation method
    translation_method = 'no_translation'

    # Change default input
    common_input.update(
        dict(branch_to_manipulate_number=branch_number, translation_method=translation_method))

    # Run with only rotation around original surface normal
    common_input.update(dict(azimuth_angle=azimuth_angle))

    # Run manipulate branch
    manipulate_branch(**common_input)

    # Set file paths
    base_path = get_path_names(common_input['input_filepath'])
    old_centerlines_path = base_path + "_centerline.vtp"
    new_centerlines_path = base_path + "_centerline_rotated.vtp"

    # Read centerlines
    old_centerlines = read_polydata(old_centerlines_path)
    old_centerlines = [extract_single_line(old_centerlines, i) for i in range(old_centerlines.GetNumberOfLines())]

    new_centerlines = read_polydata(new_centerlines_path)
    new_centerlines = [extract_single_line(new_centerlines, i) for i in range(new_centerlines.GetNumberOfLines())]

    # Compare end point of all centerlines
    unchanged_centerlines = compare_end_points(new_centerlines, old_centerlines)

    assert len(unchanged_centerlines) < len(old_centerlines)


def compare_end_points(new_centerlines, old_centerlines):
    new_unchanged_centerlines = []
    for line0 in new_centerlines:
        line0_end = get_end_point(line0)
        for line1 in old_centerlines:
            line1_end = get_end_point(line1)
            if line0_end == line1_end:
                new_unchanged_centerlines.append(line0_end)

    return new_unchanged_centerlines

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import numpy as np
import pytest

from .fixtures import surface_paths
from morphman import manipulate_bifurcation, read_command_line_bifurcation
from morphman.common.surface_operations import get_centerline_tolerance
from morphman.common.common import get_path_names, get_distance
from morphman.common.vtk_wrapper import read_polydata, get_vtk_point_locator, extract_single_line


@pytest.mark.parametrize("angle", [20 / 180 * np.pi, -20 / 180 * np.pi])
def test_bifurcation_angle(surface_paths, angle):
    common_input = read_command_line_bifurcation(surface_paths[0], surface_paths[1])
    common_input.update(dict(angle=angle,
                             region_of_interest="commandline",
                             region_points=[35.8, 59.8, 39.7, 76.8, 54.7, 53.2]))

    manipulate_bifurcation(**common_input)

    # Read in files to compute angle
    base_path = get_path_names(common_input["input_filepath"])
    old_centerlines = read_polydata(base_path + "_centerline_par.vtp")
    new_centerlines = read_polydata(base_path + "_centerline_interpolated_ang.vtp")
    end_points = read_polydata(base_path + "_clippingpoints.vtp")

    # Start points
    start_point1 = end_points.GetPoint(1)
    start_point2 = end_points.GetPoint(2)

    # Get relevant centerlines
    cl_old_1 = -1
    cl_old_2 = -1
    cl_new_1 = -1
    cl_new_2 = -1
    tol = get_centerline_tolerance(old_centerlines)

    for i in range(old_centerlines.GetNumberOfLines()):
        line_old = extract_single_line(old_centerlines, i)
        line_new = extract_single_line(new_centerlines, i)

        loc_old = get_vtk_point_locator(line_old)
        loc_new = get_vtk_point_locator(line_new)

        id1_old = loc_old.FindClosestPoint(start_point1)
        id2_old = loc_old.FindClosestPoint(start_point2)

        id1_new = loc_new.FindClosestPoint(start_point1)
        id2_new = loc_new.FindClosestPoint(start_point2)

        if get_distance(start_point1, line_old.GetPoint(id1_old)) < tol:
            cl_old_1 = i
            cl_old_id1 = id1_old
        if get_distance(start_point2, line_old.GetPoint(id2_old)) < tol:
            cl_old_2 = i
            cl_old_id2 = id2_old
        if get_distance(start_point1, line_new.GetPoint(id1_new)) < tol:
            cl_new_1 = i
            cl_new_id1 = id1_new
        if get_distance(start_point2, line_new.GetPoint(id2_new)) < tol:
            cl_new_2 = i
            cl_new_id2 = id2_new

        if -1 not in [cl_old_1, cl_old_2, cl_new_1, cl_new_2]:
            break

    # Get end points
    end_point1_old = np.array(extract_single_line(old_centerlines, cl_old_1).GetPoint(cl_old_id1 + 20))
    end_point2_old = np.array(extract_single_line(old_centerlines, cl_old_2).GetPoint(cl_old_id2 + 20))
    end_point1_new = np.array(extract_single_line(new_centerlines, cl_new_1).GetPoint(cl_new_id1 + 20))
    end_point2_new = np.array(extract_single_line(new_centerlines, cl_new_2).GetPoint(cl_new_id2 + 20))

    # Vectors
    v1_old = end_point1_old - np.array(start_point1)
    v2_old = end_point2_old - np.array(start_point2)
    v1_new = end_point1_new - np.array(start_point1)
    v2_new = end_point2_new - np.array(start_point2)

    # Normalize
    v1_old = v1_old / np.sqrt(np.sum(v1_old ** 2))
    v2_old = v2_old / np.sqrt(np.sum(v2_old ** 2))
    v1_new = v1_new / np.sqrt(np.sum(v1_new ** 2))
    v2_new = v2_new / np.sqrt(np.sum(v2_new ** 2))

    # Angle
    first_daughter_branch_angle_change = np.arccos(np.dot(v1_old, v1_new))
    second_daughter_branch_angle_change = np.arccos(np.dot(v2_old, v2_new))

    assert abs(first_daughter_branch_angle_change
               + second_daughter_branch_angle_change
               - 2 * abs(common_input["angle"])) < 0.01

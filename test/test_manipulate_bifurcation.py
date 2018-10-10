##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
from os import path
relative_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(relative_path, '..', 'src'))
sys.path.insert(0, "../src")

import pytest
import numpy as np
from .fixtures import common_input
from manipulate_bifurcation import rotate_branches
from common import get_path_names, read_polydata, get_locator, get_tolerance, \
                   extract_single_line, distance


@pytest.fixture(scope='module')
def init_data():
    # Global parameters
    l1 = l2 = False
    bif = False
    lower = True
    cylinder_factor = 7
    aneurysm = False
    anu_num = 0
    version = True

    dit = dict(basedir)
    return dic

@pytest.mark.parametrize("angle", [-20 / 180 * np.pi, 20 / 180 * np.pi])
def test_bifurcation_angle(common_input, angle):
    common_input.update(dict(keep_fixed_1 = False,
                             keep_fixed_2 = False,
                             bif = False,
                             lower = False,
                             cylinder_factor = 7,
                             angle = angle,
                             region_of_interest = "commandline",
                             region_points = [74.02971649169,
                                              54.54873275756,
                                              43.89385986328,
                                              35.75942230224,
                                              59.80244827270,
                                              39.67420196533]))

    rotate_branches(**common_input)

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
    tol = get_tolerance(old_centerlines)

    for i in range(old_centerlines.GetNumberOfLines()):
        line_old = extract_single_line(old_centerlines, i)
        line_new = extract_single_line(new_centerlines, i)

        loc_old = get_locator(line_old)
        loc_new = get_locator(line_new)

        id1_old = loc_old.FindClosestPoint(start_point1)
        id2_old = loc_old.FindClosestPoint(start_point2)

        id1_new = loc_old.FindClosestPoint(start_point1)
        id2_new = loc_old.FindClosestPoint(start_point2)

        print(distance(start_point1, old_centerlines.GetPoint(id1_old)), tol)
        if distance(start_point1, line_old.GetPoint(id1_old)) < tol:
            cl_old_1 = i
            cl_old_id1 = id1_old
        if distance(start_point2, line_old.GetPoint(id2_old)) < tol:
            cl_old_2 = i
            cl_old_id2 = id2_old
        if distance(start_point1, line_new.GetPoint(id1_new)) < tol:
            cl_new_1 = i
            cl_new_id1 = id1_new
        if distance(start_point2, line_new.GetPoint(id2_new)) < tol:
            cl_new_2 = i
            cl_new_id2 = id2_new

        if -1 not in [cl_old_1, cl_old_2, cl_new_1, cl_new_2]:
            break

    # Get end points
    end_point1_old = np.array(old_centerlines.GetPoint(cl_old_id1 + 20))
    end_point2_old = np.array(old_centerlines.GetPoint(cl_old_id2 + 20))
    end_point1_new = np.array(new_centerlines.GetPoint(cl_new_id1 + 20))
    end_point2_new = np.array(new_centerlines.GetPoint(cl_new_id2 + 20))

    # Vectors
    v1_old = end_point1_old - np.array(start_point1)
    v2_old = end_point2_old - np.array(start_point2)
    v1_new = end_point1_new - np.array(start_point1)
    v2_new = end_point2_new - np.array(start_point2)

    # Normalize
    v1_old = v1_old / np.sqrt(np.sum(v1_old**2))
    v2_old = v2_old / np.sqrt(np.sum(v2_old**2))
    v1_new = v1_new / np.sqrt(np.sum(v1_new**2))
    v2_new = v2_new / np.sqrt(np.sum(v2_new**2))

    # Angle
    angle_original = np.arccos(np.dot(v1_old, v2_old))
    angle_rotated = np.arccos(np.dot(v1_new, v2_new))

    assert angle_original - angle_rotated -common_input["angle"] < 0.1

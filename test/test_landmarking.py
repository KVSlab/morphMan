##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
from os import path
relative_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(relative_path, '..', 'src'))

import pytest
from manipulate_bend import *
from automated_geometric_quantities import compute_angle

sys.path.insert(0, '../src/')


@pytest.fixture(scope='module')
def init_data():
    # Global parameters
    basedir = "testdata"
    case = "P0134"
    point_path = "carotid_siphon_points.particles"
    dirpath = path.join(basedir, case)
    name = "surface"
    smooth_factor = 0.25
    dic = dict(point_path=point_path, dirpath=dirpath, name=name, smooth=True, smooth_factor=smooth_factor)
    return dic


def test_increase_siphon_angle(init_data):
    name = init_data['name']
    smooth = init_data['smooth']
    dirpath = init_data['dirpath']
    point_path = init_data['point_path']
    smooth_factor = init_data['smooth_factor']
    alpha = -0.1
    beta = 0.4
    method = "plane"
    move_vessel(dirpath, smooth, name, point_path, alpha, beta, smooth_factor)
    new_centerlines_path = path.join(dirpath, name, "new_centerlines_alpha_%s_beta_%s.vtp"
                                     % (alpha, beta))
    new_cl = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(dirpath, point_path, name, alpha,
                                              beta, method, new_centerline=new_cl)
    assert angle_original < angle_new

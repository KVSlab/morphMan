import pytest
from manipulate_bend import *


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

@pytest.mark.parametrize("angle", [-20, 20])
def test_bifurcation_angle(init_data, angle):
    angle = 6 # pi / 6  = 30 deg
    #TODO: Get center and end_points (before and after rotation) and measure angle

    rotate_brances(path.join(basedir, folder), name, smooth, smooth_factor, angle, l1,
                   l2, bif, lower, cylinder_factor, aneurysm, anu_num, resampling_step, version)

    angle_original = 0.5
    angle_new = 0.7
    assert angle_original < angle_new

def test_decrease_bifurcation_angle():
    angle = -6 - # pi / 6  = 30 deg
    rotate_brances(path.join(basedir, folder), name, smooth, smooth_factor, angle, l1,
                   l2, bif, lower, cylinder_factor, aneurysm, anu_num, resampling_step, version)
    angle_original = 0.5
    angle_new = 0.7
    assert angle_original > angle_new


test_increase_curvature()

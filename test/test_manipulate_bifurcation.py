import pytest
from fixture import common_input
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
def test_bifurcation_angle(common_input, angle):
    common_input.update(dict(keep_fixed_1 = False,
                             keep_fixed_2 = False,
                             bif = False,
                             lower = False,
                             cylinder_factor = 7,
                             angle = angle,
                             region_of_interest = "commandline"
                             region_points = [x, y, z, x, y, z]))

    rotate_brances(**common_input)

    # TODO: Compute new and old angles
    angle_original = 0.5
    angle_new = 0.7
    assert angle_original < angle_new

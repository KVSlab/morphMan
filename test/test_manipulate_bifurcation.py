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
    smooth = True
    basedir = "testdata"
    case = "P0134"
    name = "surface"
    l1 = l2 = False
    bif = False
    lower = True
    cylinder_factor = 7
    aneurysm = False
    anu_num = 0
    resampling_step = 0.1
    version = True
    smooth_factor = 0.25
    dirpath = path.join(basedir, case)

    dit = dict(basedir)
    return dic


def test_increase_bifurcation_angle(init_data):
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

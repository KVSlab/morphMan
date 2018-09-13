import sys; sys.path.insert(0, '../src/')

from move_siphon import *

# Global parameters
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

def test_increase_bifurcation_angle():
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

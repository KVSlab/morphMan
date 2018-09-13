import sys; sys.path.insert(0, '../src/')

from move_siphon import *

# Global parameters
smooth = True
basedir = "testdata"
case = "P0134"
name = "surface"
point_path = "carotid_siphon_points.particles"
smooth_factor = 0.25
dirpath = path.join(basedir, case)

def test_increase_curvature():
    alpha = 0.4 
    beta = -0.1
    move_vessel(dirpath, smooth, name, point_path, alpha, beta)
    curvature_original = 0.5
    curvature_new = 0.7
    assert curvature_original < curvature_new

def test_decrease_curvature():
    alpha = -0.1 
    beta = 0.4
    move_vessel(dirpath, smooth, name, point_path, alpha, beta)
    curvature_original = 0.5
    curvature_new = 0.1
    assert curvature_original > curvature_new

def test_increase_siphon_angle():
    alpha = -0.1 
    beta = 0.4
    move_vessel(dirpath, smooth, name, point_path, alpha, beta)
    angle_original = 0.5
    angle_new = 0.7
    assert angle_original < angle_new 

def test_decrease_siphon_angle():
    alpha = 0.4 
    beta = -0.1
    move_vessel(dirpath, smooth, name, point_path, alpha, beta)
    angle_original = 0.5
    angle_new = 0.7
    assert angle_original > angle_new 


test_increase_curvature()

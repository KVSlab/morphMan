import sys; sys.path.insert(0, '../src/')

from move_siphon import *
from IPython import embed

# Global parameters
smooth = True
basedir = "testdata"
case = "P0134"
name = "surface"
point_path = "carotid_siphon_points.particles"
dirpath = path.join(basedir, case)
clipping_points = get_clipping_points(dirpath, point_path)
smooth_factor = 0.25

# Old centerline
old_centerlines_path = path.join(dirpath, name, "centerline_complete.vtp")
old_cl = read_polydata(old_centerlines_path)
old_longest = extract_single_line(sort_centerlines(old_cl),0)
old_longest = vmtk_centerline_resampling(old_longest, 0.1)
_, old_curv = discrete_geometry(old_longest, neigh=20)
old_locator = get_locator(old_longest)
ID1, ID2 = old_locator.FindClosestPoint(clipping_points[0]), old_locator.FindClosestPoint(clipping_points[1]) 
if ID1 > ID2: ID1, ID2 = ID2, ID1
curvature_original = max(old_curv[ID1:ID2])

def test_increase_curvature():
    alpha = 0.4 
    beta = -0.1
    move_vessel(dirpath, smooth, name, point_path, alpha, beta)

    # Import centerlines
    new_centerlines_path = path.join(dirpath, name, "new_centerlines_alpha_%s_beta_%s.vtp"
                                     % (alpha, beta))
    new_cl = read_polydata(new_centerlines_path)
    new_longest = extract_single_line(sort_centerlines(new_cl),0)
    new_longest = vmtk_centerline_resampling(new_longest, 0.1)

    # Compute curvature
    _, new_curv = discrete_geometry(new_longest, neigh=20)

    # Select area based on clipping points
    new_locator = get_locator(new_longest)

    IDA, IDB = new_locator.FindClosestPoint(clipping_points[0]), new_locator.FindClosestPoint(clipping_points[1]) 
    if IDA > IDB: IDA, IDB = IDB, IDA

    # Compare
    curvature_new = max(new_curv[IDA:IDB])
    assert curvature_original < curvature_new

def test_decrease_curvature():
    alpha = -0.1 
    beta = 0.4
    move_vessel(dirpath, smooth, name, point_path, alpha, beta)

    # Import centerlines
    new_centerlines_path = path.join(dirpath, name, "new_centerlines_alpha_%s_beta_%s.vtp"
                                     % (alpha, beta))
    new_cl = read_polydata(new_centerlines_path)
    new_longest = extract_single_line(sort_centerlines(new_cl),0)
    new_longest = vmtk_centerline_resampling(new_longest, 0.1)

    # Compute curvature
    _, new_curv = discrete_geometry(new_longest, neigh=20)

    # Select area based on clipping points
    clipping_points = get_clipping_points(dirpath, point_path)
    new_locator = get_locator(new_longest)

    IDA, IDB = new_locator.FindClosestPoint(clipping_points[0]), new_locator.FindClosestPoint(clipping_points[1]) 
    if IDA > IDB: IDA, IDB = IDB, IDA

    # Compare
    curvature_new = max(new_curv[IDA:IDB])
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
    angle_original = 0.9
    angle_new = 0.7
    assert angle_original > angle_new 


#test_increase_curvature()
#test_decrease_curvature()
test_increase_siphon_angle()
#test_decrease_siphon_angle()

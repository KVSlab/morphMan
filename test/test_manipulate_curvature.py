import pytest
from fixtures import common_input
from manipulate_curvature import curvature_variations
from common import get_path_names, extract_single_line, read_polydata, get_locator
from automated_geometric_quantities import discrete_geometry


@pytest.mark.parameterize("smooth_line" [True, False])
def test_decrease_curvature(common_input, smooth_line):
    # Get region points
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 1)
    n = centerline.GetNumberOfPoints()
    region_points = list(centerline.GetPoint(int(n * 0.1))) + list(centerline.GetPoint(int(n * 0.4)))

    # Set problem specific parameters
    common_input.update(dict(resampling_step = 0.1,
                             smooth_factor_line = 1.5,
                             iterations = 200,
                             region_of_interest = "commandline",
                             region_points = region_points,
                             smooth_line = smooth_line))

    # Manipulate surface
    curvature_variations(**common_input)

    # Select and compare altered region
    p1 = np.asarray(region_points[:3])
    p2 = np.asarray(region_points[3:])

    old_centerlines_path = base_path + "_centerline.vtp"
    old_centerlines = read_polydata(old_centerlines_path)
    old_locator = get_locator(extract_single_line(old_centerlines, 0))
    old_id1 = old_locator.FindClosestPoint(p1)
    old_id2 = old_locator.FindClosestPoint(p2)
    old_centerline = extract_single_line(old_centerlines, 0, startID=old_id1, endID=old_id2)

    new_centerlines_path = base_path + "_centerline_new.vtp"
    new_centerlines = read_polydata(new_centerlines_path)
    new_locator = get_locator(extract_single_line(new_centerlines, 0))
    new_id1 = new_locator.FindClosestPoint(p1)
    new_id2 = new_locator.FindClosestPoint(p2)
    new_centerline = extract_single_line(new_centerlines, 0, startID=new_id1, endID=new_id2)

    # Comute curvature and assert
    _, old_curvature = discrete_geometry(old_centerline, neigh=20)
    _, new_curvature = discrete_geometry(new_centerline, neigh=20)
    old_mean_curv = np.mean(old_curvature)
    new_mean_curv = np.mean(new_curvature)
    if smooth_line:
        assert old_mean_curv > new_mean_curv
    else:
        assert old_mean_curv < new_mean_curv

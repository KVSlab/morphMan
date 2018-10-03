import pytest
from manipulate_curvature import *

relative_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(relative_path, '..', 'src'))


@pytest.fixture(scope='module')
def init_data():
    # Global parameters
    input_filepath = "./testdata/P0134/surface/model.vtp"
    output_filepath = "./testdata/P0134/surface/model_output.vtp"
    resampling_step = 0.1
    smooth_factor_line = 1.5
    iterations = 200
    smooth_factor_voro = 0.25
    poly_ball_size = [120, 120, 120]
    smooth = True
    no_smooth = False
    no_smooth_point = None

    # Get region points
    region_of_interest = "commandline"
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 1)
    n = centerline.GetNumberOfPoints()
    region_points = list(centerline.GetPoint(int(n * 0.1))) + list(centerline.GetPoint(int(n * 0.4)))

    return dict(input_filepath=input_filepath, smooth=smooth,
                smooth_factor_voro=smooth_factor_voro, smooth_factor_line=smooth_factor_line,
                iterations=iterations, output_filepath=output_filepath,
                poly_ball_size=poly_ball_size, region_of_interest=region_of_interest,
                region_points=region_points, resampling_step=resampling_step,
                no_smooth=no_smooth, no_smooth_point=no_smooth_point)


def test_decrease_curvature(init_data):
    init_data['smooth_line'] = True
    curvature_variations(**init_data)

    # Select and compare altered region
    region_points = init_data['region_points']
    p1 = np.asarray(region_points[:3])
    p2 = np.asarray(region_points[3:])
    base_path = get_path_names(init_data['input_filepath'])

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
    assert old_mean_curv > new_mean_curv


def test_increase_curvature(init_data):
    init_data['smooth_line'] = False
    curvature_variations(**init_data)

    # Select and compare altered region
    region_points = init_data['region_points']
    p1 = np.asarray(region_points[:3])
    p2 = np.asarray(region_points[3:])
    base_path = get_path_names(init_data['input_filepath'])

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
    assert old_mean_curv < new_mean_curv

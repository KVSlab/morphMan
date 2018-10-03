import pytest
from manipulate_bend import *
from automated_geometric_quantities import compute_angle

relative_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(relative_path, '..', 'src'))


@pytest.fixture(scope='module')
def init_data():
    # Global parameters
    input_filepath = "./testdata/P0134/surface/model.vtp"
    ofile = "./testdata/P0134/surface/model_output.vtp"
    region_of_interest = "landmarking"
    region_points = None
    resampling_step = 0.1
    smooth = True
    smooth_factor = 0.25
    poly_ball_size = [120, 120, 120]
    no_smooth = False
    no_smooth_point = None
    return dict(input_filepath=input_filepath, smooth=smooth,
                smooth_factor=smooth_factor,
                output_filepath=ofile, poly_ball_size=poly_ball_size,
                no_smooth=no_smooth, no_smooth_point=no_smooth_point,
                region_of_interest=region_of_interest, region_points=region_points,
                resampling_step=resampling_step)


def test_increase_siphon_angle_alpha(init_data):
    alpha, beta = -0.2, 0.0
    init_data["alpha"] = alpha
    init_data["beta"] = beta
    input_filepath = init_data["input_filepath"]

    # Perform manipulation
    move_vessel(**init_data)

    # Compute angle
    method = "plane"
    new_centerlines_path = input_filepath[:-4] + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(input_filepath, alpha,
                                              beta, method, new_centerlines=new_centerlines)
    assert angle_original < angle_new


def test_decrease_siphon_angle_alpha(init_data):
    alpha, beta = 0.4, 0.0
    init_data["alpha"] = alpha
    init_data["beta"] = beta
    input_filepath = init_data["input_filepath"]

    # Perform manipulation
    move_vessel(**init_data)

    # Compute angle
    method = "plane"
    new_centerlines_path = input_filepath[:-4] + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(input_filepath, alpha,
                                              beta, method, new_centerlines=new_centerlines)
    assert angle_original > angle_new


def test_increase_siphon_angle_beta(init_data):
    alpha, beta = 0.0, 0.4
    init_data["alpha"] = alpha
    init_data["beta"] = beta
    input_filepath = init_data["input_filepath"]

    # Perform manipulation
    move_vessel(**init_data)

    # Compute angle
    method = "plane"
    new_centerlines_path = input_filepath[:-4] + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(input_filepath, alpha,
                                              beta, method, new_centerlines=new_centerlines)
    assert angle_original < angle_new


def test_decrease_siphon_angle_beta(init_data):
    alpha, beta = 0.0, - 0.15
    init_data["alpha"] = alpha
    init_data["beta"] = beta
    input_filepath = init_data["input_filepath"]

    # Perform manipulation
    move_vessel(**init_data)

    # Compute angle
    method = "plane"
    new_centerlines_path = input_filepath[:-4] + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(input_filepath, alpha,
                                              beta, method, new_centerlines=new_centerlines)
    assert angle_original > angle_new

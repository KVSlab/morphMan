import pytest
import sys
from os import path

rel_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(rel_path, '..', 'src'))
from manipulate_bend import *
from automated_geometric_quantities import compute_angle


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


def test_increase_siphon_angle(init_data):
    alpha, beta = 0.0, 0.4
    init_data["alpha"] = alpha
    init_data["beta"] = beta
    input_filepath = init_data["input_filepath"]
    # move_vessel(input_filepath, smooth, smooth_factor, alpha, beta, output_filepath, poly_ball_size, no_smooth,
    #            no_smooth_point)

    move_vessel(**init_data)

    method = "plane"
    new_centerlines_path = input_filepath[:-4] + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_cl = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(input_filepath, alpha,
                                              beta, method, new_centerlines=new_cl)
    assert angle_original < angle_new

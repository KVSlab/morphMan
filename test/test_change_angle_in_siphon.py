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
    base_path = "../testdata/P0134/surface/model.vtp"
    smooth = True
    smooth_factor = 0.25
    output_filepath = "../testdata/P0134/surface/model_output.vtp"
    poly_ball_size = [120, 120, 120]
    no_smooth = False
    no_smooth_point = None
    init_data_dict = dict(base_path=base_path, smooth_factor=smooth_factor, smooth=smooth,
                          output_filepath=output_filepath, poly_ball_size=poly_ball_size,
                          no_smooth_point=no_smooth_point, no_smooth=no_smooth)
    return init_data_dict


def test_increase_siphon_angle(init_data):
    smooth = init_data['smooth']
    input_filepath = init_data['base_path']
    smooth_factor = init_data['smooth_factor']
    no_smooth = init_data['no_smooth']
    no_smooth_point = init_data['no_smooth_point']
    output_filepath = init_data['output_filepath']
    alpha = 0.0
    beta = 0.4
    poly_ball_size = init_data['poly_ball_size']

    move_vessel(input_filepath, smooth, smooth_factor, alpha, beta, output_filepath, poly_ball_size, no_smooth,
                no_smooth_point)

    method = "plane"
    new_centerlines_path = input_filepath[:-4] + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_cl = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(input_filepath, alpha,
                                              beta, method, new_centerline=new_cl)
    assert angle_original < angle_new

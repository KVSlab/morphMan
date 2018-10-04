import pytest
from manipulate_bend import move_vessel
from automated_geometric_quantities import compute_angle
from fixtures import common_input
from common import read_polydata, get_path_names

@pytest.mark.parametrize("alpha,beta",
                         [(-0.2,  0.0),
                          ( 0.2,  0.0),
                          ( 0.0,  0.2),
                          ( 0.0, -0.2),
                          ( 0.2, -0.2)])
def test_siphon_angle(common_input, alpha, beta):
    # Set problem specific parameters
    common_input.update(dict(alpha=alpha,
                             beta=beta,
                             region_of_interest="commandline",
                             # TODO: Set points
                             region_points=[x, y, z, x, y, z],
                             resampling_step = 0.1))

    # Perform manipulation
    move_vessel(**common_input)

    # Compute angle
    base_path = get_path_names(common_input['input_filepath'])
    new_centerlines_path = base_path + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(input_filepath, alpha, beta, "plane",
                                              new_centerlines)

    if alpha < 0 and beta == 0 or beta > 0 and alpha == 0:
        assert angle_original < angle_new
    elif alpha > 0 and beta == 0 or beta < 0 and alpha == 0:
        assert angle_original > angle_new
    else:
        assert angle_original > angle_new


@pytest.mark.parametrize("alpha,beta",
                         [(-0.2,  0.0),
                          ( 0.4,  0.0),
                          ( 0.0,  0.4),
                          ( 0.0, -0.2),
                          ( 0.2, -0.2)])
def test_siphon_curvature(common_input, alpha, beta):
    # Set problem specific parameters
    common_input.update(dict(alpha=alpha,
                             beta=beta,
                             region_of_interest="commandline",
                             # TODO: Set points
                             region_points=[x, y, z, x, y, z],
                             resampling_step = 0.1))

    # Perform manipulation
    move_vessel(**common_input)

    # Compute angle
    base_path = get_path_names(common_input['input_filepath'])
    new_centerlines_path = base_path + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_curvature(input_filepath, alpha, beta, "disc",
                                                  new_centerlines)

    if alpha < 0 and beta == 0 or beta < 0 and alpha == 0:
        assert angle_original < angle_new
    elif alpha > 0 and beta == 0 or beta > 0 and alpha == 0:
        assert angle_original > angle_new
    else:
        assert angle_original > angle_new

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import sys
from os import path
relative_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(relative_path, '..', 'src'))

import pytest
from manipulate_bend import move_vessel
from automated_geometric_quantities import compute_angle, compute_curvature
from fixtures import common_input
from common import read_polydata, get_path_names

@pytest.mark.parametrize("alpha,beta",
                         [(-0.2,  0.0),
                          ( 0.2,  0.0),
                          ( 0.0,  0.2),
                          ( 0.0, -0.2),
                          ( 0.2, -0.2)])
def test_siphon(common_input, alpha, beta):
    # Set problem specific parameters
    common_input.update(dict(alpha=alpha,
                             beta=beta,
                             region_of_interest="commandline",
                             # TODO: Set points
                             region_points=[44.17085266113281,
                                            38.514854431152344,
                                            41.20818328857422,
                                            43.242130279541016,
                                            42.68572235107422,
                                            38.65191650390625],
                             resampling_step = 0.1))

    # Perform manipulation
    print(common_input.items())
    move_vessel(**common_input)

    # Compute angle
    base_path = get_path_names(common_input['input_filepath'])
    new_centerlines_path = base_path + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(common_input["input_filepath"], alpha, beta,
                                              "plane", new_centerlines,
                                              region_points=common_input["region_points"])

    if alpha < 0 and beta == 0 or beta > 0 and alpha == 0:
        assert angle_original < angle_new
    elif alpha > 0 and beta == 0 or beta < 0 and alpha == 0:
        assert angle_original > angle_new
    else:
        assert angle_original > angle_new

    # Compute angle
    #angle_new, angle_original = compute_curvature(common_input["input_filepath"], alpha,
    #                                              beta, "disc", new_centerlines,
    #                                              region_points=common_input["region_points"])

    #if alpha < 0 and beta == 0 or beta < 0 and alpha == 0:
    #    assert angle_original < angle_new
    #elif alpha > 0 and beta == 0 or beta > 0 and alpha == 0:
    #    assert angle_original > angle_new
    #else:
        assert angle_original > angle_new

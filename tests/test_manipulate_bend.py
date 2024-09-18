##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

import pytest

from morphman import manipulate_bend, read_command_line_bend
from morphman.common.tools_common import get_path_names
from morphman.common.vtk_wrapper import read_polydata
from morphman.misc import compute_angle

from .fixtures import surface_paths


@pytest.mark.parametrize(
    "alpha,beta", [(-0.2, 0.0), (0.2, 0.0), (0.0, 0.2), (0.0, -0.2), (0.2, -0.2)]
)
def test_siphon(surface_paths, alpha, beta):
    # Set problem specific parameters
    common_input = read_command_line_bend(surface_paths[0], surface_paths[1])
    common_input.update(
        dict(
            alpha=alpha,
            beta=beta,
            region_of_interest="commandline",
            region_points=[
                44.17085266113281,
                38.514854431152344,
                41.20818328857422,
                43.242130279541016,
                42.68572235107422,
                38.65191650390625,
            ],
        )
    )

    # Perform manipulation
    manipulate_bend(**common_input)

    # Compute angle
    base_path = get_path_names(common_input["input_filepath"])
    new_centerlines_path = base_path + "_centerlines_alpha_%s_beta_%s.vtp" % (
        alpha,
        beta,
    )
    new_centerlines = read_polydata(new_centerlines_path)
    angle_new, angle_original = compute_angle(
        common_input["input_filepath"],
        alpha,
        beta,
        "plane",
        new_centerlines,
        region_of_interest=common_input["region_of_interest"],
        region_points=common_input["region_points"],
    )

    if alpha < 0 and beta == 0 or beta > 0 and alpha == 0:
        assert angle_original < angle_new
    elif alpha > 0 and beta == 0 or beta < 0 and alpha == 0:
        assert angle_original > angle_new
    else:
        assert angle_original > angle_new

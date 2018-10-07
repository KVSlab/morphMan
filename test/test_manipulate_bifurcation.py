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
import numpy as np
from fixtures import common_input
from manipulate_bifurcation import rotate_branches


@pytest.fixture(scope='module')
def init_data():
    # Global parameters
    l1 = l2 = False
    bif = False
    lower = True
    cylinder_factor = 7
    aneurysm = False
    anu_num = 0
    version = True

    dit = dict(basedir)
    return dic

@pytest.mark.parametrize("angle", [-20 / 180 * np.pi, 20 / 180 * np.pi])
def test_bifurcation_angle(common_input, angle):
    common_input.update(dict(keep_fixed_1 = False,
                             keep_fixed_2 = False,
                             bif = False,
                             lower = False,
                             cylinder_factor = 7,
                             angle = angle,
                             region_of_interest = "commandline",
                             region_points = [74.02971649169,
                                              54.54873275756,
                                              43.89385986328,
                                              35.75942230224,
                                              59.80244827270,
                                              39.67420196533]))

    #if angle < 0:
    #    common_input["input_filepath"] = common_input["input_filepath"].replace(".vtp",
    #                                                                            "2.vtp")
    rotate_branches(**common_input)

    # TODO: Compute new and old angles
    angle_original = 0.5
    angle_new = 0.7
    assert angle_original < angle_new

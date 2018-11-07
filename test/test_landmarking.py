##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.


import sys
from os import path
relative_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(relative_path, '..', 'morphman'))
sys.path.insert(0, "../morphman")

import pytest
from common import *

from automated_landmarking import landmarking_bogunovic, landmarking_piccinelli
from .fixtures import common_input

print(relative_path)


@pytest.mark.parametrize("algorithm", ["bogunovic", "piccinelli"])
def test_landmarking(common_input, algorithm):
    # Get region points
    base_path = get_path_names(common_input["input_filepath"])
    relevant_outlets = [35.8, 59.8, 39.7, 76.8, 54.7, 53.2]
    ica_centerline = extract_ica_centerline(base_path, common_input["resampling_step"],
                                            relevant_outlets=relevant_outlets)

    landmark_input = dict(
        centerline=ica_centerline, base_path=base_path, curv_method="spline", algorithm=algorithm,
        resampling_step=common_input["resampling_step"], smooth_line=False, nknots=8, iterations=100)

    if algorithm == "bogunovic":
        landmark_input.update(smoothing_factor=1.0)
        landmarks = landmarking_bogunovic(**landmark_input)
        assert len(landmarks) == 4
    elif algorithm == "piccinelli":
        landmark_input.update(smoothing_factor_torsion=1.0, smoothing_factor_curv=1.0)
        landmarks = landmarking_piccinelli(**landmark_input)
        assert len(landmarks) > 0


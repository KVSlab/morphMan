##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.


import pytest

from .fixtures import surface_paths
from morphman.common.common import get_path_names
from morphman.common.surface_operations import extract_ica_centerline
from morphman.automated_landmarking import landmarking_bogunovic, landmarking_piccinelli


@pytest.mark.parametrize("algorithm", ["bogunovic", "piccinelli"])
def test_automated_landmarking(surface_paths, algorithm):
    # Get region points
    base_path = get_path_names(surface_paths[0])
    relevant_outlets = [35.8, 59.8, 39.7, 76.8, 54.7, 53.2]
    ica_centerline = extract_ica_centerline(base_path, surface_paths[0],  0.1, relevant_outlets=relevant_outlets)

    landmark_input = dict(
        centerline=ica_centerline, base_path=base_path, approximation_method="spline", algorithm=algorithm,
        resampling_step=0.1, smooth_line=False, nknots=8, iterations=100)

    if algorithm == "bogunovic":
        landmark_input.update(smoothing_factor=1.0, coronal_axis='z')
        landmarks = landmarking_bogunovic(**landmark_input)
        assert len(landmarks) == 4

    if algorithm == "piccinelli":
        landmark_input.update(smoothing_factor_torsion=1.0, smoothing_factor_curv=1.0)
        landmarks = landmarking_piccinelli(**landmark_input)
        assert len(landmarks) > 0

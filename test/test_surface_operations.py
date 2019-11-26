##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

# noinspection PyUnresolvedReferences
from IPython import embed

from .fixtures import surface_paths
from morphman.common.surface_operations import *


def test_relevant_outlets(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    surface = read_polydata(input_filepath)

    # Provide fake outlets
    params = {'relevant_outlet_1': [0.1, 0.1, 0.1], 'relevant_outlet_2': [1, 1, 1]}
    write_parameters(params, base_path)

    outlets = get_relevant_outlets(surface, base_path)

    assert len(outlets) == 2

    assert outlets[0] == [0.1, 0.1, 0.1]
    assert outlets[1] == [1, 1, 1]


def test_compute_centers(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    surface = read_polydata(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    inlet, outlets = compute_centers(surface, base_path)

    assert len(inlet) // 3 == 1
    assert len(outlets) // 3 == centerlines.GetNumberOfLines()



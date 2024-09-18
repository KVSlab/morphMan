##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.
## This software is distributed WITHOUT ANY WARRANTY; without even
## the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE. See the above copyright notices for more information.

import sys

import pytest

from morphman.common.surface_operations import *
from morphman.common.tools_common import get_path_names
from morphman.common.vtk_wrapper import read_polydata

from .fixtures import surface_paths


@pytest.mark.skipif(sys.platform == "win32", reason="does not run on windows")
def test_uncapp_surface(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)
    capped = vmtk_cap_polydata(surface)

    uncapped = get_uncapped_surface(capped)

    is_capped, _ = is_surface_capped(uncapped)
    assert not is_capped


def test_relevant_outlets(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    surface = read_polydata(input_filepath)

    # Provide fake outlets
    params = {"relevant_outlet_1": [0.1, 0.1, 0.1], "relevant_outlet_2": [1, 1, 1]}
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


def test_compute_area_ratio(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    area_ratio, center = compute_circleness(surface)

    assert type(area_ratio) == np.float64
    assert len(center) // 3 == 1


def test_capped_surface(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)
    capped = vmtk_cap_polydata(surface)

    is_capped, number_outlets = is_surface_capped(surface)
    assert not is_capped

    is_capped, _ = is_surface_capped(capped)
    assert is_capped

    assert number_outlets == 7


def test_if_surface_is_merged(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    surface = read_polydata(input_filepath)
    capped = vmtk_cap_polydata(surface)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    # Sanity check on input model
    check_if_surface_is_merged(capped, centerlines, input_filepath)


def test_compute_voronoi_diagram(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    surface = read_polydata(input_filepath)
    capped = vmtk_cap_polydata(surface)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    voronoi = prepare_voronoi_diagram(
        capped, centerlines, base_path, False, 0, True, None, None, None, 0.1
    )

    assert voronoi.GetNumberOfPoints() > 0


def test_compute_centerlines(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    surface = read_polydata(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    inlet, outlets = compute_centers(surface, base_path)

    assert len(inlet) // 3 == 1
    assert len(outlets) // 3 == centerlines.GetNumberOfLines()


def test_prepare_surface(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)

    surface, capped_surface = prepare_surface(base_path, input_filepath)

    is_closed = is_surface_capped(capped_surface)

    assert is_closed

    assert surface.GetNumberOfCells() > 0

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from .fixtures import surface_paths
from morphman.common.surface_operations import *


def test_remove_distant_points(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")
    voronoi = read_polydata(base_path + "_voronoi.vtp")

    voronoi_removed = remove_distant_voronoi_points(voronoi, centerlines)

    assert voronoi_removed.GetNumberOfPoints() < voronoi.GetNumberOfPoints()


def test_smooth_voronoi_diagram(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")
    voronoi = read_polydata(base_path + "_voronoi.vtp")

    voronoi_smoothed = smooth_voronoi_diagram(voronoi, centerlines, 0.95)

    assert voronoi_smoothed.GetNumberOfPoints() == voronoi.GetNumberOfPoints()


def test_create_new_surface(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    surface = read_polydata(input_filepath)
    voronoi = read_polydata(base_path + "_voronoi.vtp")

    new_surface = create_new_surface(voronoi, poly_ball_size=[75, 75, 75])

    # Sanity check - Compare high resolution surface with low poly version
    assert new_surface.GetNumberOfCells() > 0
    assert new_surface.GetNumberOfCells() < surface.GetNumberOfCells()


def test_split_voronoi_diagram(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    voronoi = read_polydata(base_path + "_voronoi.vtp")

    centerlines = read_polydata(base_path + "_centerline.vtp")

    centerline1 = extract_single_line(centerlines, 0)
    centerline2 = extract_single_line(centerlines, 1)

    split_voronoi = get_split_voronoi_diagram(voronoi, [centerline1, centerline2])

    assert len(split_voronoi) == 2

    assert split_voronoi[0].GetNumberOfCells() > split_voronoi[1].GetNumberOfCells()

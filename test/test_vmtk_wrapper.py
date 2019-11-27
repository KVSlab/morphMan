##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.
from IPython import embed

from .fixtures import surface_paths
from morphman.common.surface_operations import *
from morphman.common.vmtk_wrapper import *
import os


def test_compute_geometric_features(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    centerline = vmtk_compute_geometric_features(centerline, False, outputsmoothed=False, factor=1.0, iterations=100)

    features = get_number_of_arrays(centerline)[1]

    assert len(features) == 6
    assert "Curvature" in features
    assert "Torsion" in features
    assert radiusArrayName in features
    assert "FrenetTangent" in features
    assert "FrenetNormal" in features
    assert "FernetBiNormal" in features


def test_smooth_centerline(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    smoothed_centerline = vmtk_smooth_centerline(centerline, 100, 1.0)

    centerline = vmtk_compute_geometric_features(centerline, False)
    smoothed_centerline = vmtk_compute_geometric_features(smoothed_centerline, False)

    curvature = get_point_data_array("Curvature", centerline)
    curvature_smooth = get_point_data_array("Curvature", smoothed_centerline)

    assert np.mean(curvature_smooth) < np.mean(curvature)


def test_resample_centerline(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    centerline_resampled = vmtk_resample_centerline(centerline, 1.0)
    centerline_resampled_more = vmtk_resample_centerline(centerline, 5.0)

    assert centerline_resampled.GetNumberOfPoints() < centerline.GetNumberOfPoints()

    assert centerline_resampled_more.GetNumberOfPoints() < centerline_resampled.GetNumberOfPoints()


def test_compute_sections(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    centerline = vmtk_resample_centerline(centerline, 1.0)
    surface = read_polydata(input_filepath)

    _, sections = vmtk_compute_centerline_sections(surface, centerline)

    assert sections.GetNumberOfCells() > 0


def test_compute_attributes(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)

    centerline = vmtk_compute_centerline_attributes(centerline)

    attributes = get_number_of_arrays(centerline)[1]

    assert parallelTransportNormalsArrayName in attributes
    assert abscissasArrayName in attributes
    assert radiusArrayName in attributes


def test_cap_polydata(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    capped_surface = vmtk_cap_polydata(surface)

    assert capped_surface.GetNumberOfCells() > surface.GetNumberOfCells()


def test_smooth_surface(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    smooth_surface = vmtk_smooth_surface(surface, "laplace")

    assert smooth_surface.GetNumberOfCells() == surface.GetNumberOfCells()


def test_compute_voronoi_diagram(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)
    voronoi_path = "voronoi_test.vtp"

    voronoi = vmtk_compute_voronoi_diagram(surface, voronoi_path)

    assert os.path.exists(voronoi_path)

    assert voronoi.GetNumberOfPoints() > 0

    os.remove(voronoi_path)


def test_surface_connectivity(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    connected_surface = vmtk_surface_connectivity(surface).Surface

    # Sanity check on connected model
    assert connected_surface.GetNumberOfCells() == surface.GetNumberOfCells()


def test_endpoint_extractor(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    extractor = vmtk_endpoint_extractor(centerlines, 0)
    clipped_centerlines = extractor.Centerlines

    assert clipped_centerlines.GetNumberOfLines() > centerlines.GetNumberOfLines()


def test_compute_surface_normals(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    surface_normals = vmtk_compute_surface_normals(surface)
    properties = get_number_of_arrays(surface_normals)[1]

    assert surfaceNormalsArrayName in properties


def test_branch_extractor(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    centerlines_branched = vmtk_compute_branch_extractor(centerlines)

    assert centerlines.GetNumberOfLines() < centerlines_branched.GetNumberOfLines()


def test_compute_surface_distance(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    distance_surface = vmtk_surface_distance(surface, surface)

    # Sanity test with same surface
    properties = get_number_of_arrays(distance_surface)[1]

    assert "Distance" in properties

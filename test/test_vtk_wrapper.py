##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.
from IPython import embed

from .fixtures import surface_paths
from morphman.common.surface_operations import *
from morphman.common.vtk_wrapper import *
from morphman.common.vmtk_wrapper import *
import os


def test_read_polydata(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    assert centerlines is not None


def test_write_polydata(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")
    save_path = "centerlines_test.vtp"

    write_polydata(centerlines, save_path)

    assert os.path.exists(save_path)

    os.remove(save_path)


def test_number_of_arrays(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    # Get MISR array
    n, names = get_number_of_arrays(centerlines)

    assert n == 1
    assert radiusArrayName in names


def test_extract_single_line(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    # Get MISR array
    first_line = extract_single_line(centerlines, 0)

    assert first_line.GetNumberOfLines() == 1
    assert first_line.GetNumberOfPoints() < centerlines.GetNumberOfPoints()


def test_merge_polydata(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")
    n_lines = centerlines.GetNumberOfLines()

    double_centerlines = vtk_merge_polydata([centerlines, centerlines])

    assert double_centerlines.GetNumberOfLines() == 2 * n_lines


def test_compute_normals(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    cell_normals = vtk_compute_polydata_normals(surface, compute_cell_normals=True)

    assert cell_normals.GetNumberOfCells() > 0


def test_compute_gradients(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    cell_normals = vtk_compute_polydata_normals(surface, compute_cell_normals=True)
    gradients = vtk_compute_normal_gradients(cell_normals)

    assert gradients.GetNumberOfCells() > 0


def test_get_cell_data_array(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    cell_normals = vtk_compute_polydata_normals(surface, compute_cell_normals=True)
    gradients = vtk_compute_normal_gradients(cell_normals)

    gradient_array = get_cell_data_array("Gradients", gradients, k=9)

    assert type(gradient_array) == np.ndarray

    assert len(gradient_array) == gradients.GetNumberOfCells()


def test_get_point_data_array(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    misr = get_point_data_array(radiusArrayName, centerlines)

    assert len(misr) == centerlines.GetNumberOfPoints()


def test_clean_surface(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    clean_surface = vtk_clean_polydata(surface)

    assert clean_surface.GetNumberOfCells() == surface.GetNumberOfCells()


def test_surface_connectivity(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    # Sanity check using already connected surface
    connected_surface = vtk_compute_connectivity(surface, mode="Largest")

    assert connected_surface.GetNumberOfCells() == surface.GetNumberOfCells()


def test_vtk_threshold(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    cell_normals = vtk_compute_polydata_normals(surface, compute_cell_normals=True)
    gradients = vtk_compute_normal_gradients(cell_normals)

    gradients_threshold = vtk_compute_threshold(gradients, "Gradients")

    assert gradients_threshold.GetNumberOfCells() < gradients.GetNumberOfCells()


def test_get_boundary_edges(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    edges = vtk_extract_feature_edges(surface)

    assert edges.GetNumberOfCells() > 0


def test_get_vtk_array(surface_paths):
    size = 50
    vtk_array = get_vtk_array(radiusArrayName, 1, size)

    assert vtk_array.GetSize() == size
    assert vtk_array.GetName() == radiusArrayName


def test_create_vtk_array(surface_paths):
    values = np.linspace(0, 50)
    size = len(values)
    vtk_array = create_vtk_array(values, radiusArrayName)

    assert vtk_array.GetSize() == size
    assert vtk_array.GetName() == radiusArrayName


def test_vtk_point_locator(surface_paths):
    input_filepath = surface_paths[0]
    base_path = get_path_names(input_filepath)
    centerlines = read_polydata(base_path + "_centerline.vtp")

    point_locator = get_vtk_point_locator(centerlines)

    p0 = centerlines.GetPoint(0)

    closest_point_id = point_locator.FindClosestPoint(p0)

    p = centerlines.GetPoint(closest_point_id)

    assert p == p0


def test_triangulate_surface(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    triangulated_surface = vtk_triangulate_surface(surface)

    assert triangulated_surface.GetNumberOfCells() == surface.GetNumberOfCells()


def test_compute_area(surface_paths):
    input_filepath = surface_paths[0]
    surface = read_polydata(input_filepath)

    mass = vtk_compute_mass_properties(surface)

    assert mass > 0


def test_compute_plane():
    origin = [2, 3, -1]
    n = [3.14, -5, 1]

    plane = vtk_plane(origin, n)

    assert list(plane.GetNormal()) == n
    assert list(plane.GetOrigin()) == origin


def test_get_sphere():
    center = [1, 2, 1]
    radius = 5.5

    sphere = vtk_sphere(center, radius)

    assert sphere.GetRadius() == radius
    assert list(sphere.GetCenter()) == center

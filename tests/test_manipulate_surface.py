##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.
## This software is distributed WITHOUT ANY WARRANTY; without even
## the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE. See the above copyright notices for more information.

import numpy as np

from morphman import manipulate_surface, read_command_line_surface
from morphman.common.vmtk_wrapper import vmtk_surface_curvature
from morphman.common.voronoi_operations import create_new_surface
from morphman.common.vtk_wrapper import (
    get_point_data_array,
    read_polydata,
    write_polydata,
)

from .fixtures import surface_paths


def test_smooth_surface(surface_paths):
    # Get default inputs
    common_input = read_command_line_surface(surface_paths[0], surface_paths[1])
    common_input["poly_ball_size"] = [250, 250, 250]

    # Manipulate surface
    manipulate_surface(**common_input)

    old_voronoi = read_polydata(surface_paths[0].replace(".vtp", "_voronoi.vtp"))
    old_surface = create_new_surface(
        old_voronoi, poly_ball_size=common_input["poly_ball_size"]
    )
    new_surface = read_polydata(surface_paths[1])

    # Compare curvature
    old_surface = vmtk_surface_curvature(
        old_surface,
        curvature_type="gaussian",
        absolute=True,
        median_filtering=True,
        bounded_reciporcal=True,
    )
    new_surface = vmtk_surface_curvature(
        new_surface,
        curvature_type="gaussian",
        absolute=True,
        median_filtering=True,
        bounded_reciporcal=True,
    )
    write_polydata(new_surface, "tmp_new.vtp")
    write_polydata(old_surface, "tmp_old.vtp")
    old_curvature = get_point_data_array("Curvature", old_surface, k=1)
    new_curvature = get_point_data_array("Curvature", new_surface, k=1)

    assert np.mean(new_curvature) < np.mean(old_curvature)


def test_noise_surface(surface_paths):
    # Get default inputs
    common_input = read_command_line_surface(surface_paths[0], surface_paths[1])
    common_input.update(
        dict(
            noise=True,
            smooth=False,
            frequency=5,
            frequency_deviation=0.01,
            add_noise_lower_limit=0.8,
            add_noise_upper_limit=0.95,
            radius_min=1.1,
            radius_max=1.35,
            poly_ball_size=[250, 250, 250],
        )
    )

    # Manipulate surface
    manipulate_surface(**common_input)

    old_surface = read_polydata(surface_paths[0])
    new_surface = read_polydata(surface_paths[1])

    # Compare curvature
    old_surface = vmtk_surface_curvature(
        old_surface,
        curvature_type="gaussian",
        absolute=True,
        median_filtering=True,
        bounded_reciporcal=True,
    )
    new_surface = vmtk_surface_curvature(
        new_surface,
        curvature_type="gaussian",
        absolute=True,
        median_filtering=True,
        bounded_reciporcal=True,
    )

    old_curvature = get_point_data_array("Curvature", old_surface, k=1)
    new_curvature = get_point_data_array("Curvature", new_surface, k=1)

    assert np.mean(new_curvature) > np.mean(old_curvature)

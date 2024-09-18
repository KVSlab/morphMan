## Writen by Aslak W. Bergersen, 2019
## Modified by Henrik A. Kjeldsberg, 2022
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from morphman.common import *


def landmark_atrium(input_path, resampling_step):
    print("-- Extracting left atrium and the left atrial appendage")
    extract_left_atrium_and_appendage(input_path, resampling_step)

    print("-- Separating the left atrium and the left atrial appendage")
    separate_left_atrium_and_appendage(input_path, resampling_step)

    print("-- Landmarking of the left atrium complete")


def extract_left_atrium_and_appendage(input_path, resampling_step):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along enterlines.

    Args:
        input_path (str): Path to folder with model files
        resampling_step (float): Resampling length for centerline resampling.
    """
    print("--- Load model file\n")
    clipped_model = input_path.replace(".stl", "_la_and_laa.vtp")

    new_cl = input_path.replace(".stl", "_centerline.vtp")
    surface = read_polydata(input_path)

    capped_surface = vmtk_cap_polydata(surface)

    # Centers
    inlet, outlets = compute_centers(surface)
    p_outlet = np.array(inlet)

    # Centerline
    capped_surface = vtk_clean_polydata(capped_surface)
    la_centerlines, _, _ = compute_centerlines(
        inlet, outlets, new_cl, capped_surface, resampling=resampling_step, smooth=True
    )

    # Clip PVs
    N_PVs = 4
    for i in range(N_PVs):
        print(f"--- Clipping PV ({i + 1})")
        la_centerline_i = extract_single_line(la_centerlines, i)
        start = int(la_centerline_i.GetNumberOfPoints() * 0.5)
        stop = int(la_centerline_i.GetNumberOfPoints() * 0.95)

        line = extract_single_line(la_centerline_i, 0, start_id=start, end_id=stop)
        la_l = get_curvilinear_coordinate(line)
        step = 5 * np.mean(la_l[1:] - la_l[:-1])
        line = vmtk_resample_centerline(line, step)
        line = compute_splined_centerline(line, nknots=10, isline=True)
        area, sections = vmtk_compute_centerline_sections(capped_surface, line)

        # Get arrays
        a = get_point_data_array("CenterlineSectionArea", area)
        n = get_point_data_array("FrenetTangent", area, k=3)

        # Compute 'derivative' of the area
        dAdX = np.gradient(a.T[0], step)

        # Find the largest change in cross-sectional area
        lim = -50
        tol = 2
        stop_id = np.nonzero(dAdX < lim)[0][-1] + tol
        normal = n[stop_id]
        center = area.GetPoint(stop_id)

        # Clip the model
        plane = vtk_plane(center, normal)
        surface, clipped = vtk_clip_polydata(surface, plane)

        # Find part to keep
        surface = vtk_clean_polydata(surface)
        clipped = vtk_clean_polydata(clipped)
        p_boundary = np.array(
            la_centerline_i.GetPoint(la_centerline_i.GetNumberOfPoints() - 1)
        )

        surf_loc = get_vtk_point_locator(surface)
        clip_loc = get_vtk_point_locator(clipped)
        id_surf = surf_loc.FindClosestPoint(p_boundary)
        id_clip = clip_loc.FindClosestPoint(p_boundary)
        p_surface = np.array(surface.GetPoint(id_surf))
        p_clipped = np.array(clipped.GetPoint(id_clip))
        dist_surface = np.linalg.norm(p_surface - p_boundary)
        dist_clipped = np.linalg.norm(p_clipped - p_boundary)
        if dist_surface < dist_clipped:
            surface, clipped = clipped, surface

        surface = attach_clipped_regions_to_surface(surface, clipped, center)

    # Clip MV
    print("--- Clipping MV")
    # Compute 'derivative' of the area
    la_centerline_0 = extract_single_line(la_centerlines, 0)
    start = 0
    stop = int(la_centerline_0.GetNumberOfPoints() * 0.5)

    line = extract_single_line(la_centerline_0, 0, start_id=start, end_id=stop)
    la_l = get_curvilinear_coordinate(line)
    step = 5 * np.mean(la_l[1:] - la_l[:-1])
    line = vmtk_resample_centerline(line, step)
    line = compute_splined_centerline(line, nknots=10, isline=True)
    area, sections = vmtk_compute_centerline_sections(capped_surface, line)

    # Get arrays
    a = get_point_data_array("CenterlineSectionArea", area)
    n = get_point_data_array("FrenetTangent", area, k=3)

    dAdX = np.gradient(a.T[0], step)

    # Find the largest change in cross-sectional area
    lim = 50
    tol = 5
    stop_id = np.nonzero(dAdX > lim)[0][0] + tol
    normal = n[stop_id]
    center = area.GetPoint(stop_id)

    # Clip the model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)

    # Find part to keep
    surface = vtk_clean_polydata(surface)
    clipped = vtk_clean_polydata(clipped)
    p_boundary = p_outlet

    surf_loc = get_vtk_point_locator(surface)
    clip_loc = get_vtk_point_locator(clipped)
    id_surf = surf_loc.FindClosestPoint(p_boundary)
    id_clip = clip_loc.FindClosestPoint(p_boundary)
    p_surface = np.array(surface.GetPoint(id_surf))
    p_clipped = np.array(clipped.GetPoint(id_clip))
    dist_surface = np.linalg.norm(p_surface - p_boundary)
    dist_clipped = np.linalg.norm(p_clipped - p_boundary)

    if dist_surface < dist_clipped:
        surface, clipped = clipped, surface

    print("--- Saving LA and LAA to: {}".format(clipped_model))
    surface = attach_clipped_regions_to_surface(surface, clipped, center)
    write_polydata(surface, clipped_model)


def separate_left_atrium_and_appendage(input_path, resampling_step):
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along centerlines.

    Args:
        input_path (str): Path to the surface for landmarking
        resampling_step (float): Resampling length for centerline resampling.
    """
    print("--- Load model file\n")

    laa_model_path = input_path.replace(".vtp", "_laa.vtp")
    laa_centerline_path = input_path.replace(".vtp", "_laa_centerline.vtp")
    la_model_path = input_path.replace(".vtp", "_la.vtp")

    surface = read_polydata(input_path)

    print("-- Reading centerlines")
    capped_surface = vmtk_cap_polydata(surface)

    # Centers
    inlet, outlets = compute_centers(surface)

    # Pick LAA point
    p_laa = provide_point(surface)

    laa_centerline, _, _ = compute_centerlines(
        inlet,
        p_laa,
        laa_centerline_path,
        capped_surface,
        resampling=resampling_step,
        smooth=True,
    )
    id_start = int(laa_centerline.GetNumberOfPoints() * 0.1)
    id_stop = int(laa_centerline.GetNumberOfPoints() * 0.8)

    line = extract_single_line(laa_centerline, 0, start_id=id_start, end_id=id_stop)
    laa_l = get_curvilinear_coordinate(line)
    step = np.mean(laa_l[1:] - laa_l[:-1])
    line = vmtk_resample_centerline(line, step)
    line = compute_splined_centerline(line, nknots=10, isline=True)
    area, sections = vmtk_compute_centerline_sections(capped_surface, line)

    # Get arrays
    print("-- Extract point arrays")
    a = get_point_data_array("CenterlineSectionArea", area)
    n = get_point_data_array("FrenetTangent", area, k=3)

    # Stopping criteria
    dAdX = np.gradient(a.T[0], step)
    lim = -200
    tol = 5
    ID = -1

    stop_id = np.nonzero(dAdX < lim)[0][ID] + tol
    normal = n[stop_id]
    center = area.GetPoint(stop_id)

    # Clip the model
    plane = vtk_plane(center, normal)
    surface, clipped = vtk_clip_polydata(surface, plane)

    # Find part to keep
    surface = vtk_clean_polydata(surface)
    clipped = vtk_clean_polydata(clipped)
    p_boundary = p_laa

    surf_loc = get_vtk_point_locator(surface)
    clip_loc = get_vtk_point_locator(clipped)
    id_surf = surf_loc.FindClosestPoint(p_boundary)
    id_clip = clip_loc.FindClosestPoint(p_boundary)
    p_surface = np.array(surface.GetPoint(id_surf))
    p_clipped = np.array(clipped.GetPoint(id_clip))
    dist_surface = np.linalg.norm(p_surface - p_boundary)
    dist_clipped = np.linalg.norm(p_clipped - p_boundary)

    # Extract LAA:
    if dist_surface > dist_clipped:
        surface, clipped = clipped, surface

    surface_whole = attach_clipped_regions_to_surface(surface, clipped, center)
    laa_surface = get_surface_closest_to_point(surface_whole, center)

    # Extract LA:
    if dist_surface < dist_clipped:
        surface, clipped = clipped, surface

    surface_whole = attach_clipped_regions_to_surface(surface, clipped, center)
    la_surface = get_surface_closest_to_point(surface_whole, center)
    print("--- Saving LAA to: {}".format(laa_model_path))
    write_polydata(laa_surface, laa_model_path)

    print("--- Saving LA to: {}".format(la_model_path))
    write_polydata(la_surface, la_model_path)


def get_surface_closest_to_point(clipped, point):
    """Check the connectivty of a clipped surface, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        clipped (vtkPolyData): The clipped segments of the surface.
        point (list): The point of interest. Keep region closest to this point

    Returns:
        surface (vtkPolyData): The surface where only one segment has been removed.
    """
    connectivity = vtk_compute_connectivity(clipped, mode="All")
    if connectivity.GetNumberOfPoints() == 0:
        return clipped

    region_id = get_point_data_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(
            vtk_compute_threshold(
                connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0
            )
        )
        locator = get_vtk_point_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(point))
        distances.append(get_distance(region_point, point))

    # Remove the region with the closest distance
    region_of_interest = regions[distances.index(min(distances))]

    return region_of_interest


def read_command_line_landmark(input_path=None):
    """
    THe landmarking performs two tasks: 'extract' for extracting the left atrium and appendage,
    and 'separate' for separating the left atrium and appendage.

    Args:
        input_path (str): Input file path, positional argument with default None.
    """
    description = "Landmarks the left atrium and appendage of a 3D model."
    parser = ArgumentParser(
        description=description, formatter_class=RawDescriptionHelpFormatter
    )

    required = input_path is None
    add_common_arguments(parser, required=required)

    if required:
        args = parser.parse_args()
    else:
        args = parser.parse_args(["-i", input_path])

    return dict(input_path=args.input_path, resampling_step=args.resampling_step)


def main_landmark():
    landmark_atrium(**read_command_line_landmark())


if __name__ == "__main__":
    landmark_atrium(**read_command_line_landmark())

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

import math
from os import makedirs

from vtk.util import numpy_support

from morphman.common.voronoi_operations import *


def get_relevant_outlets(surface, base_path):
    """
    Extract relevant outlets of the
    input surface model.

    Args:
        surface(vtkPolyData): Surface model.
        base_path (str): Location of info-file.

    Returns:
        relevant_outlets (list): List of relevant outlet IDs.
    """
    # Check if info exists
    if not path.isfile(base_path + "_info.json"):
        provide_relevant_outlets(surface, base_path)

    # Open info
    parameters = get_parameters(base_path)
    relevant_outlets = []
    for key, value in list(parameters.items()):
        if key.startswith("relevant_outlet_"):
            relevant_outlets.append(value)

    if not relevant_outlets:
        relevant_outlets = provide_relevant_outlets(surface, base_path)

    return relevant_outlets


def compute_centers(polydata, case_path=None):
    """
    Compute the center of all the openings in the surface. The inlet is chosen based on
    the largest area.

    Args:
        polydata (vtkPolyData): centers of the openings
        case_path (str): path to case directory.

    Returns:
        inlet (list): A list of points.
        outlet (list): A flattened list with all the outlets.
    """
    # Get cells which are open
    cells = vtk_extract_feature_edges(polydata)

    if cells.GetNumberOfCells() == 0:
        print("WARNING: The model is capped, so it is uncapped, but the method is experimental.")
        uncapped_surface = get_uncapped_surface(polydata)
        compute_centers(uncapped_surface, case_path)

    # Compute connectivity of the cells
    outputs = vtk_compute_connectivity(cells)

    # Get connectivity array
    region_array = get_point_data_array("RegionId", outputs)

    # Get points
    points = np.zeros((region_array.shape[0], 3))
    for i in range(region_array.shape[0]):
        points[i,] = outputs.GetPoint(i)

    # Get area and center
    area = []
    center = []
    for i in range(int(region_array.max()) + 1):
        # Extract points for this opening
        tmp_points = points[region_array[:, 0] == i, :]

        # Extract the surface edge
        boundary_points = vtk_compute_threshold(outputs, "RegionId", lower=i - 0.1, upper=i + 0.1,
                                                threshold_type="between", source=0)

        # Create surface for computing area of opening
        delaunay = vtk.vtkDelaunay2D()
        delaunay.SetInputData(boundary_points)
        delaunay.Update()

        # Add quanteties
        area.append(vtk_compute_mass_properties(delaunay.GetOutput()))
        center.append(np.mean(tmp_points, axis=0))

    # Store the center and area
    inlet_ind = area.index(max(area))
    if case_path is not None:
        info = {"inlet": center[inlet_ind].tolist(), "inlet_area": area[inlet_ind]}
        p = 0
        for i in range(len(area)):
            if i == inlet_ind:
                p = -1
                continue

            info["outlet%d" % (i + p)] = center[i].tolist()
            info["outlet%s_area" % (i + p)] = area[i]

        write_parameters(info, case_path)

    inlet_center = center[inlet_ind].tolist()
    center.pop(inlet_ind)

    center_ = [item for sublist in center for item in sublist]

    return inlet_center, center_


def provide_relevant_outlets(surface, dir_path=None):
    """
    Get relevant outlets from user
    selected points on a input surface.

    Args:
        surface (vtkPolyData): Surface model.
        dir_path (str): Location of info.json file

    Returns:
        points (list): List of relevant outlet IDs
    """

    # Fix surface
    cleaned_surface = vtk_clean_polydata(surface)
    triangulated_surface = vtk_triangulate_surface(cleaned_surface)

    # Select seeds
    print("-- Please select the two relevant outlets in the interactive window.")
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(triangulated_surface)
    seed_selector.text = "Please select the two relevant outlets, \'u\' to undo\n"
    seed_selector.Execute()

    point_seed_ids = seed_selector.GetTargetSeedIds()
    get_point = surface.GetPoints().GetPoint
    points = [list(get_point(point_seed_ids.GetId(i))) for i in range(point_seed_ids.GetNumberOfIds())]
    info = {}

    if dir_path is not None:
        for i in range(len(points)):
            info["relevant_outlet_%d" % i] = points[i]
        write_parameters(info, dir_path)

    return points


def get_inlet_and_outlet_centers(surface, base_path, flowext=False):
    """Get the centers of the inlet and outlets.

    Args:
        surface (vtkPolyData): An open surface.
        base_path (str): Path to the case file.
        flowext (bool): Turn on/off flow extension.

    Returns:
        inlet (list): A flatt list with the point of the inlet
        outlet (list): A flatt list with the points of all the outlets.
    """
    # Check if info exists
    if flowext or not path.isfile(base_path + "_info.json"):
        compute_centers(surface, base_path)

    # Open info
    parameters = get_parameters(base_path)
    outlets = []
    inlet = []
    for key, value in list(parameters.items()):
        if key == "inlet":
            inlet = value
        elif "outlet" in key and "area" not in key and "relevant" not in key:
            outlets += value

    num_outlets = len(outlets) // 3
    if num_outlets != 0:
        outlets = []
        for i in range(num_outlets):
            outlets += parameters["outlet%d" % i]

    if inlet == [] and outlets == []:
        inlet, outlets = compute_centers(surface, base_path)

    return inlet, outlets


def get_clipped_capped_surface(surface, centerlines, clipspheres=0):
    """A method for clipping a capped outlets. The branches will be clipped some distance
    from the outlets.

    Args:
        surface (vtkPolyData): Surface to clipp
        centerlines (vtkPolyData): Centerlines to mark the in and outlets.
        clipspheres (float): Number of end point spheres

    Returns:
        surface (vtkPolyData): Clipped surface
    """
    extractor = vmtk_endpoint_extractor(centerlines, clipspheres)
    clipped_centerlines = extractor.Centerlines

    clipper = vmtk_branch_clipper(clipped_centerlines, surface)
    surface = clipper.Surface

    connector = vmtk_surface_connectivity(surface)
    surface = connector.Surface

    return surface


def compute_circleness(surface):
    """Compute the area ratio betwen minimum circle and the maximum circle.

    Args:
        surface (vtkPolyData): Boundary edges of an opening

    Returns:
        circleness (float): Area ratio
        center (list): Center of the opening.
    """
    edges = vtk_extract_feature_edges(surface)

    # Get points
    points = []
    for i in range(edges.GetNumberOfPoints()):
        points.append(edges.GetPoint(i))

    # Compute center
    points = np.array(points)
    center = np.mean(np.array(points), axis=0)

    # Compute ratio between max inscribed sphere, and min inscribed "area"
    point_radius = np.sqrt(np.sum((points - center) ** 2, axis=1))
    argsort = np.argsort(point_radius)
    if point_radius[argsort[1]] / point_radius[argsort[0]] > 5:
        radius_min = point_radius[argsort[1]]
    else:
        radius_min = point_radius.min()

    min_area = math.pi * radius_min ** 2
    max_area = math.pi * point_radius.max() ** 2
    circleness = max_area / min_area

    return circleness, center


def is_surface_capped(surface):
    """Checks if the surface is closed, and how many openings there are.

    Args:
        surface (vtkPolyData): Surface to be checked

    Returns:
        open (boolean): Open or closed surface
        number (int): Number of integer
    """
    # Get boundary cells
    cells = vtk_extract_feature_edges(surface)
    if cells.GetNumberOfCells() == 0:
        return True, 0
    else:
        outlets = vtk_compute_connectivity(cells, mode="All")
        number = get_point_data_array("RegionId", outlets).max()
        return number == 0, int(number)


def get_uncapped_surface(surface, gradients_limit=0.15, area_limit=0.3, circleness_limit=3):
    """
    A rule-based method for removing endcapps on a surface. The method considers the
    gradient of the normals, the size of the region, and how similar it is to a circle.

    Args:
        surface (vtkPolyData): Surface to be uncapped.
        gradients_limit (float): Upper limit for gradients of normals.
        area_limit (float): Lower limit of the area.
        circleness_limit (float): Upper limit of the circleness.

    Returns:
        surface (vtkPolyData): The uncapped surface.

    """

    cell_normals = vtk_compute_polydata_normals(surface, compute_cell_normals=True)

    gradients = vtk_compute_normal_gradients(cell_normals)

    # Compute the magnitude of the gradient
    gradients_array = get_cell_data_array("Gradients", gradients, 9)
    gradients_magnitude = np.sqrt(np.sum(gradients_array ** 2, axis=1))

    # Mark all cells with a gradient magnitude less then gradient_limit
    end_capp_array = gradients_magnitude < gradients_limit
    end_capp_vtk = get_vtk_array("Gradients_mag", 1, end_capp_array.shape[0])
    for i, p in enumerate(end_capp_array):
        end_capp_vtk.SetTuple(i, [p])
    gradients.GetCellData().AddArray(end_capp_vtk)

    # Extract capps
    end_capps = vtk_compute_threshold(gradients, "Gradients_mag", lower=0.5, upper=1.5,
                                      threshold_type="between", source=1)

    # Get connectivity
    end_capps_connectivity = vtk_compute_connectivity(end_capps)
    region_array = get_point_data_array("RegionId", end_capps_connectivity)

    # Compute area for each region
    area = []
    circleness = []
    regions = []
    centers_edge = []
    limit = 0.1
    for i in range(int(region_array.max()) + 1):
        regions.append(vtk_compute_threshold(end_capps_connectivity, "RegionId", lower=(i - limit),
                                             upper=(i + limit), threshold_type="between", source=0))
        circ, center = compute_circleness(regions[-1])
        circleness.append(circ)
        centers_edge.append(center)
        area.append(vtk_compute_mass_properties(regions[-1]))

    # Only keep outlets with circleness < circleness_limit and area > area_limit
    circleness_ids = np.where(np.array(circleness) < circleness_limit)
    region_ids = np.where(np.array(area) > area_limit)
    regions = [regions[i] for i in region_ids[0] if i in circleness_ids[0]]
    centers_edge = [centers_edge[i] for i in region_ids[0] if i in circleness_ids[0]]

    # Mark the outlets on the original surface
    mark_outlets = create_vtk_array(np.zeros(surface.GetNumberOfCells()), "outlets", k=1)
    locator = get_vtk_cell_locator(surface)
    tmp_center = [0, 0, 0]
    for region in regions:
        centers_filter = vtk.vtkCellCenters()
        centers_filter.SetInputData(region)
        centers_filter.VertexCellsOn()
        centers_filter.Update()
        centers = centers_filter.GetOutput()

        for i in range(centers.GetNumberOfPoints()):
            centers.GetPoint(i, tmp_center)
            p = [0, 0, 0]
            cell_id = vtk.mutable(0)
            sub_id = vtk.mutable(0)
            dist = vtk.mutable(0)
            locator.FindClosestPoint(tmp_center, p, cell_id, sub_id, dist)
            mark_outlets.SetTuple(cell_id, [1])

    surface.GetCellData().AddArray(mark_outlets)

    # Remove the outlets from the original surface
    uncapped_surface = vtk_compute_threshold(surface, "outlets", lower=0, upper=0.5, threshold_type="between", source=1)

    # Check if some cells where not marked
    remove = True
    while remove:
        locator = get_vtk_cell_locator(uncapped_surface)
        mark_outlets = create_vtk_array(np.zeros(uncapped_surface.GetNumberOfCells()), "outlets", k=1)
        remove = False
        for center in centers_edge:
            locator.FindClosestPoint(center, p, cell_id, sub_id, dist)
            if dist < 0.01:
                remove = True
                mark_outlets.SetTuple(cell_id, [1])

        uncapped_surface.GetCellData().AddArray(mark_outlets)

        if remove:
            uncapped_surface = vtk_compute_threshold(uncapped_surface, "outlets", lower=0,
                                                     upper=0.5, threshold_type="between", source=1)

    return uncapped_surface


def check_if_surface_is_merged(surface, centerlines, output_filepath):
    """
    Check if surface has overlapping regions.

    Args:
        surface (vtkPolyData): Surface model.
        centerlines (vtkPolyData): New centerlines.
        output_filepath (str): Filepath of output model.
    """
    # Check if the manipulated centerline and the centerline from the new surface
    # significantly differ, if so it is likely that part of the surface is now merged
    centerlines = vmtk_resample_centerline(centerlines, length=0.1)
    inlet = centerlines.GetPoint(0)
    outlets = []
    lines_to_compare = []
    for i in range(centerlines.GetNumberOfLines()):
        lines_to_compare.append(extract_single_line(centerlines, i))
        outlets += lines_to_compare[-1].GetPoint(lines_to_compare[-1].GetNumberOfPoints() - 1)
    lines_to_check, _, _ = compute_centerlines(inlet, outlets, None, surface,
                                               resampling=0.1, recompute=True)

    for i in range(centerlines.GetNumberOfLines()):
        line_to_compare = vmtk_resample_centerline(lines_to_compare[i], length=0.1)
        line_to_check = vmtk_resample_centerline(extract_single_line(lines_to_check, i), length=0.1)
        # Compare distance between points along both centerliens
        n = min([line_to_check.GetNumberOfPoints(), line_to_compare.GetNumberOfPoints()])
        tolerance = get_centerline_tolerance(line_to_compare) * 500
        for j in range(n):
            p1 = np.asarray(line_to_check.GetPoint(j))
            p2 = np.asarray(line_to_compare.GetPoint(j))
            dist = get_distance(p1, p2)
            if dist > tolerance:
                tmp_path = output_filepath.replace(".vtp", "_ERROR_MERGED.vtp")
                write_polydata(surface, tmp_path)
                raise RuntimeError(("\nERROR: Model has most likely overlapping regions." +
                                    " Please check the surface model {} and provide other" +
                                    " parameters for the manipulation or" +
                                    " poly_ball_size.").format(tmp_path))


def prepare_output_surface(surface, original_surface, new_centerline, output_filepath,
                           test_merge=False, changed=False, old_centerline=None,
                           removed=[[1e9, 1e9, 1e9]]):
    """After manipulation preparing the surface for output. This method clipps the
    outlets, slightly smooths the surface, and (potentially) tests if the surface is is
    merged.

    Args:
        surface (vtkPolyData): The new surface after manipulation.
        original_surface (vtkPolyData): The original surface inputed for manipulation.
        new_centerline (vtkPolyData): The centerline after manipulation.
        output_filepath (str): The user-defined path to the output.
        test_merge (bool): Turn on/off testing if the surface is merged.
        changed (bool): If the manipulated surface has changed the location of the
        inlet/outlet.
        old_centerline (vtkPolyData): The old centerline for the original centerline.

    Returns:
        surface (vtkPolyData): The surface ready for output.
    """
    # Check if the folder for the output exits
    if not path.exists(path.dirname(output_filepath)):
        if path.dirname(output_filepath) != "":
            makedirs(path.dirname(output_filepath))

    # Get planes if outlets of the original surface
    boundary_edges = vtk_extract_feature_edges(original_surface)
    boundary_connectivity = vtk_compute_connectivity(boundary_edges)

    vtk_array = boundary_connectivity.GetPointData().GetArray("RegionId")
    vtk_points = boundary_connectivity.GetPoints().GetData()
    region_id = numpy_support.vtk_to_numpy(vtk_array)
    points = numpy_support.vtk_to_numpy(vtk_points)

    centerline = new_centerline if old_centerline is None else old_centerline
    outlets = []
    lines = []
    for i in range(centerline.GetNumberOfLines()):
        lines.append(extract_single_line(centerline, i))
        outlets.append(lines[-1].GetPoint(lines[-1].GetNumberOfPoints() - 1))
    inlet_point = lines[-1].GetPoint(0)

    if changed and old_centerline is None:
        print("WARNING: The changed flag is true, but the old centerline is not provided," +
              " and the outlet location can therefore not be changed.")

    # Get information from the original geometry
    inlet = False
    for i in range(region_id.max() + 1):
        # Get relevant points
        tmp_points = points[region_id == i]

        # Get normal
        tmp_normal = np.cross(tmp_points[0] - tmp_points[-1],
                              tmp_points[0] - tmp_points[tmp_points.shape[0] // 2])
        normal = tmp_normal / np.sqrt(np.sum(tmp_normal ** 2))

        # Get Center
        center = np.mean(tmp_points, axis=0)

        # Check if branch has been removed
        if np.sqrt(np.sum(np.array(removed) - center) ** 2) < 0.5:
            continue

        # Get corresponding centerline to in/outlet
        if np.sqrt(np.sum((np.array(inlet_point) - center) ** 2)) < 0.5:
            line = lines[0]
            line_id = 0
            inlet = True
        else:
            line_id = np.argmin(np.sqrt(np.sum((np.array(outlets) - center) ** 2, axis=1)))
            line = lines[line_id]

        # Set correct direction of normal
        if inlet:
            in_dir = np.array(line.GetPoint(5)) - \
                     np.array(line.GetPoint(0))
        else:
            in_dir = np.array(line.GetPoint(line.GetNumberOfPoints() - 5)) - \
                     np.array(line.GetPoint(line.GetNumberOfPoints() - 1))

        in_dir = in_dir / np.sqrt(np.sum(in_dir ** 2))
        angle = np.arccos(np.dot(in_dir, normal)) * 180 / np.pi
        normal = -normal if 90 < angle < 270 else normal

        # Mapp the old center and normals to the altered model
        if changed and old_centerline is not None:
            new_line = extract_single_line(new_centerline, line_id)

            # Set correct direction of normal
            if inlet:
                new_outlet = np.array(new_line.GetPoint(0))
                in_dir_new = np.array(new_line.GetPoint(5)) - new_outlet
                translation = new_outlet - np.array(inlet_point)
            else:
                new_outlet = np.array(new_line.GetPoint(new_line.GetNumberOfPoints() - 1))
                in_dir_new = np.array(new_line.GetPoint(new_line.GetNumberOfPoints() - 5)) - new_outlet
                translation = new_outlet - np.array(outlets[line_id])

            center += translation
            in_dir_new = in_dir_new / np.sqrt(np.sum(in_dir_new ** 2))
            in_dir_normal = np.cross(in_dir_new, in_dir)
            dir_angle = np.arccos(np.dot(in_dir, in_dir_new)) * 180 / np.pi

            translation = vtk.vtkTransform()
            translation.RotateWXYZ(-dir_angle, in_dir_normal)
            tmp_normal = normal
            normal = [0, 0, 0]
            translation.TransformNormal(tmp_normal, normal)

        # Set plane
        plane = vtk_plane(center, normal)

        # Clip data (naivly)
        surface, clipped = vtk_clip_polydata(surface, plane)

        # Reattach data which should not have been clipped
        surface = attach_clipped_regions_to_surface(surface, clipped, center)
        inlet = False

    # Perform a 'light' smoothing to obtain a nicer surface
    surface = vmtk_smooth_surface(surface, method="laplace", iterations=100)

    # Clean surface
    surface = vtk_clean_polydata(surface)
    surface = vtk_triangulate_surface(surface)

    # Capped surface
    capped_surface = vmtk_cap_polydata(surface)
    if test_merge:
        check_if_surface_is_merged(capped_surface, new_centerline, output_filepath)

    return surface


def attach_clipped_regions_to_surface(surface, clipped, center):
    """Check the connectivty of a clipped surface, and attach all sections which are not
    closest to the center of the clipping plane.

    Args:
        surface (vtkPolyData):
        clipped (vtkPolyData): The clipped segments of the surface.
        center (list): The center of the clipping point

    Returns:
        surface (vtkPolyData): The surface where only one segment has been removed.
    """
    connectivity = vtk_compute_connectivity(clipped, mode="All")
    if connectivity.GetNumberOfPoints() == 0:
        return surface
    region_id = get_point_data_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(vtk_compute_threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0))
        locator = get_vtk_point_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(center))
        distances.append(get_distance(region_point, center))

    # Remove the region with the closest distance
    regions.pop(distances.index(min(distances)))

    # Add the other regions back to the surface
    surface = vtk_merge_polydata(regions + [surface])
    surface = vtk_clean_polydata(surface)
    surface = vtk_triangulate_surface(surface)

    return surface


def prepare_voronoi_diagram(capped_surface, centerlines, base_path, smooth, smooth_factor, no_smooth, no_smooth_point,
                            voronoi, pole_ids, resampling_length, absolute=False, upper=None):
    """
    Compute and smooth voronoi diagram of surface model.

    Args:
        capped_surface (polydata): Cappedsurface model to create a Voronoi diagram of.
        base_path (str): Absolute path to surface model path.
        voronoi (vtkPolyData): Voronoi diagram.
        pole_ids (vtkIDList): Pole ids of Voronoi diagram.
        smooth (bool): Voronoi is smoothed if True.
        smooth_factor (float): Smoothing factor for voronoi smoothing.
        centerlines (vtkPolyData): Centerlines throughout geometry.
        no_smooth (bool): Part of Voronoi is not smoothed.
        no_smooth_point (vtkPolyData): Point which defines unsmoothed area.
        resampling_length (float): Length of resampling the centerline.
        absolute (bool): Turn on/off absolute values for the smoothing. Default is off.
        upper (int): Set an upper limit for the smoothing factor. Default is None.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram of surface.
    """
    # Check if a region should not be smoothed
    if smooth and no_smooth:
        no_smooth_cl = get_no_smooth_cl(capped_surface, centerlines, base_path, smooth, no_smooth, voronoi,
                                        no_smooth_point, pole_ids, resampling_length)
    else:
        no_smooth_cl = None

    if voronoi is None:
        voronoi = vmtk_compute_voronoi_diagram(capped_surface, base_path + "_voronoi.vtp")

    # Smooth voronoi
    voronoi_smoothed_path = base_path + "_voronoi_smoothed.vtp"
    surface_smoothed_path = base_path + "_smoothed.vtp"
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi = smooth_voronoi_diagram(voronoi, centerlines, smooth_factor, no_smooth_cl)
        write_polydata(voronoi, voronoi_smoothed_path)

        # Create new surface from the smoothed Voronoi
        surface_smoothed = create_new_surface(voronoi)
        write_polydata(surface_smoothed, surface_smoothed_path)
    elif smooth:
        voronoi = read_polydata(voronoi_smoothed_path)

    return voronoi


def compute_centerlines(inlet, outlet, filepath, surface, resampling=1.0, smooth=False,
                        num_iter=100, smooth_factor=0.1, end_point=1, method="pointlist",
                        recompute=False, voronoi=None, pole_ids=None, base_path=None):
    """Wrapper for vmtkcenterlines and vmtkcenterlinesmoothing.

    Args:
        inlet (list): point of the inlet
        outlet (list): flatt list of the outlet points
        filepath (str): path to where to store the centerline
        surface (vtkPolyData): surface to get the centerline from.
        resampling (float): resampling step length.
        smooth (bool): smooth centerline or not.
        num_iter (int): number of iterations in smooth.
        smooth_factor (float): smoothing factor.
        end_point (int): 0 or 1, include end point in centerline.
        method (str): method for setting the inlet and outlet location
        recompute (bool): if filepath exists, but the centerline should be computed again
        anyway.
        voronoi (vtkPolyData): Optional argument for setting the Voronoi diagram.
        pole_ids (vtkIdList): A vtkIdList coupling the surface with the voronoi diagram
        base_path (str): path to the case

    Returns:
        centerline (vtkPolyData): centerline of the surface.
        voronoi (vtkPolyData): Voronoi data.
        pole_ids (vtkIdList): vtkIdList coupling the surface and the voronoi diagram.
    """
    if path.isfile(str(filepath)) and not recompute:  # Filepath might be None
        if base_path is not None and path.isfile(base_path + "_voronoi.vtp"):
            voronoi = read_polydata(base_path + "_voronoi.vtp")
            pole_ids = read_polydata(base_path + "_pole_ids.np", datatype="vtkIdList")
        else:
            voronoi = None
            pole_ids = None

        return read_polydata(filepath), voronoi, pole_ids

    centerlines, centerlines_output = vmtk_compute_centerlines(end_point, inlet, method, outlet, pole_ids, resampling,
                                                               surface, voronoi)

    if smooth:
        centerlines_output = vmtk_smooth_centerline(centerlines_output, num_iter, smooth_factor)

    # Save the computed centerline.
    if filepath is not None:
        write_polydata(centerlines_output, filepath)

    voronoi = centerlines.VoronoiDiagram
    pole_ids = centerlines.PoleIds
    if base_path is not None:
        write_polydata(voronoi, base_path + "_voronoi.vtp")
        write_polydata(pole_ids, base_path + "_pole_ids.np", datatype="vtkIdList")

    return centerlines_output, voronoi, pole_ids


def prepare_surface(base_path, surface_path):
    """
    Clean and check connectivity of surface.
    Capps or uncapps surface at inlet and outlets.

    Args:
        base_path (str): Absolute path to base folder.
        surface_path (str): Path to surface.

    Returns:
        open_surface (vtkPolyData): Open surface.
    Returns:
        capped_surface (vtkPolyData): Closed surface.
    """
    # Check if surface path exists
    surface_capped_path = base_path + "_capped.vtp"
    if not path.exists(surface_path):
        RuntimeError("Could not find the file: {}".format(surface_path))

    # Clean surface
    surface = read_polydata(surface_path)
    surface = vtk_clean_polydata(surface)
    surface = vtk_triangulate_surface(surface)

    # Check connectivity and only choose the surface with the largest area
    parameters = get_parameters(base_path)
    if "check_surface" not in parameters.keys():
        connected_surface = vtk_compute_connectivity(surface, mode="Largest")
        if connected_surface.GetNumberOfCells() != surface.GetNumberOfCells():
            write_polydata(surface, surface_path.replace(".vtp", "_unconnected.vtp"))
            write_polydata(connected_surface, surface_path)
            surface = connected_surface

        parameters["check_surface"] = True
        write_parameters(parameters, base_path)

    # Get a capped and uncapped version of the surface
    cap_bool, num_out = is_surface_capped(surface)
    if cap_bool:
        open_surface = get_uncapped_surface(surface)
        cap_bool, num_out = is_surface_capped(open_surface)
        print(("WARNING: Tried to automagically uncapp the input surface. Uncapped {}" +
               " inlet/outlets in total. If this number if incorrect please provide an" +
               " uncapped surface as input, use the clipp_capped_surface" +
               " method, or vmtksurfaceendclipper.").format(num_out))
        capped_surface = surface
        write_polydata(capped_surface, surface_capped_path)
        write_polydata(open_surface, surface_path)
    else:
        open_surface = surface
        if path.exists(surface_capped_path):
            capped_surface = read_polydata(surface_capped_path)
        else:
            capped_surface = vmtk_cap_polydata(surface)
            write_polydata(capped_surface, surface_capped_path)

    return open_surface, capped_surface


def extract_ica_centerline(base_path, input_filepath, resampling_step, relevant_outlets=None):
    """
    Extract a centerline from the inlet to the first branch.

    Args:
        base_path (str): Path to the case folder.
        input_filepath (str): Path to the original model.
        resampling_step (float): Resampling step length of the extracted centerline.
        relevant_outlets (ndarray): Array containing points corresponding to outlets

    Returns:
        centerline (vtkPolyData): Extracted centerline.
    """
    # TODO: Extract ICA centerline by comparing daughter branch cross-section areas
    centerlines_path = base_path + "_centerline.vtp"
    ica_centerline_path = base_path + "_ica.vtp"
    centerline_relevant_outlets_path = base_path + "_centerline_relevant_outlets_landmark.vtp"

    if path.exists(ica_centerline_path):
        return read_polydata(ica_centerline_path)

    # Prepare surface and identify in/outlets
    surface, capped_surface = prepare_surface(base_path, input_filepath)
    inlet, outlets = get_inlet_and_outlet_centers(surface, base_path)

    if relevant_outlets is not None:
        outlet1, outlet2 = relevant_outlets[:3], relevant_outlets[3:]
        surface_locator = get_vtk_point_locator(capped_surface)
        id1 = surface_locator.FindClosestPoint(outlet1)
        id2 = surface_locator.FindClosestPoint(outlet2)
        outlet1 = capped_surface.GetPoint(id1)
        outlet2 = capped_surface.GetPoint(id2)
    else:
        outlet1, outlet2 = get_relevant_outlets(capped_surface, base_path)
    outlets, outlet1, outlet2 = get_sorted_outlets(outlets, outlet1, outlet2, base_path)

    # Compute / import centerlines
    compute_centerlines(inlet, outlets, centerlines_path, capped_surface, resampling=resampling_step, smooth=False,
                        base_path=base_path)

    # Get relevant centerlines
    centerline_relevant_outlets = compute_centerlines(inlet, outlet1 + outlet2, centerline_relevant_outlets_path,
                                                      capped_surface, resampling=resampling_step)

    # Extract ICA centerline
    tmp_line_1 = extract_single_line(centerline_relevant_outlets[0], 0)
    tmp_line_2 = extract_single_line(centerline_relevant_outlets[0], 1)
    tolerance = get_centerline_tolerance(tmp_line_1)
    line = extract_single_line(tmp_line_1, 0, start_id=0,
                               end_id=get_diverging_point_id(tmp_line_1, tmp_line_2, tolerance))
    write_polydata(line, ica_centerline_path)

    return line


def get_no_smooth_cl(capped_surface, centerlines, base_path, smooth, no_smooth, voronoi,
                     no_smooth_point, pole_ids, resampling_length, region_points=None):
    """
    Extract a section where the Voronoi should not be smoothed
     Args:
        capped_surface (polydata): Cappedsurface model to create a Voronoi diagram of.
        centerlines (vtkPolyData): Centerlines throughout geometry.
        base_path (str): Absolute path to surface model path.
        smooth (bool): Voronoi is smoothed if True.
        no_smooth (bool): Part of Voronoi is not smoothed.
        voronoi (vtkPolyData): Voronoi diagram.
        no_smooth_point (vtkPolyData): Point which defines unsmoothed area.
        pole_ids (vtkIDList): Pole ids of Voronoi diagram.
        resampling_length (float): Length of resampling the centerline.
        region_points (list): Flatten list with the region points
     Returns:
        no_smooth_cl (vtkPolyData): Centerline section where the Voronoi should not be smoothed
    """

    no_smooth_path = base_path + "_centerline_no_smooth.vtp"
    # Get inlet and outlets
    tol = get_centerline_tolerance(centerlines)
    inlet = extract_single_line(centerlines, 0)
    inlet = inlet.GetPoint(0)
    outlets = []
    parameters = get_parameters(base_path)

    if "no_smooth_point_1" in parameters.keys():
        counter = 1
        while "no_smooth_point_{}".format(counter) in parameters.keys():
            outlets += parameters["no_smooth_point_{}".format(counter)]
            counter += 1

    elif no_smooth_point is None:
        seed_selector = vmtkPickPointSeedSelector()
        seed_selector.SetSurface(capped_surface)
        seed_selector.text = "Please place a point on the segments you do not want" + \
                             " to smooth, e.g. an aneurysm, \'u\' to undo\n"
        seed_selector.Execute()
        point_ids = seed_selector.GetTargetSeedIds()
        for i in range(point_ids.GetNumberOfIds()):
            parameters["no_smooth_point_{}".format(i)] = capped_surface.GetPoint(point_ids.GetId(i))
            outlets += capped_surface.GetPoint(point_ids.GetId(i))

    else:
        locator = get_vtk_point_locator(capped_surface)
        for i in range(len(no_smooth_point) // 3):
            tmp_id = locator.FindClosestPoint(no_smooth_point[3 * i:3 * (i + 1)])
            parameters["no_smooth_point_{}".format(i)] = capped_surface.GetPoint(tmp_id)
            outlets += capped_surface.GetPoint(tmp_id)

    # Store parameters
    write_parameters(parameters, base_path)

    # Create the centerline
    no_smooth_centerlines, _, _ = compute_centerlines(inlet, outlets, None,
                                                      capped_surface,
                                                      resampling=resampling_length,
                                                      smooth=False, voronoi=voronoi, pole_ids=pole_ids)

    # Remove the no_smooth_centerline outside the region of interest
    if region_points is not None:
        no_smooth_centerlines_list = []
        tol = get_centerline_tolerance(no_smooth_centerlines)

        for i in range(no_smooth_centerlines.GetNumberOfLines()):
            line = extract_single_line(no_smooth_centerlines, i)
            locator = get_vtk_point_locator(line)
            id1 = locator.FindClosestPoint(region_points[:3])
            id2 = locator.FindClosestPoint(region_points[3:])

            distance1 = get_distance(region_points[:3], line.GetPoint(id1))
        distance2 = get_distance(region_points[3:], line.GetPoint(id2))

        if distance1 < tol * 5 != distance2 < tol * 5:
            no_smooth_centerlines_list.append(extract_single_line(line, i, start_id=id1))

        no_smooth_centerlines = vtk_merge_polydata(no_smooth_centerlines_list)

    # Extract the centerline region which diverges from the existing centerlines
    no_smooth_segments = []
    for i in range(no_smooth_centerlines.GetNumberOfLines()):
        tmp_line = extract_single_line(no_smooth_centerlines, i)
        div_ids = []
        for j in range(centerlines.GetNumberOfLines()):
            div_ids.append(get_diverging_point_id(tmp_line, extract_single_line(centerlines, j), tol))
        div_id = max(div_ids)
        no_smooth_segments.append(extract_single_line(tmp_line, 0, start_id=div_id))

    no_smooth_cl = vtk_merge_polydata(no_smooth_segments)
    write_polydata(no_smooth_cl, no_smooth_path)
    # else:
    #    no_smooth_cl = read_polydata(no_smooth_path)

    return no_smooth_cl

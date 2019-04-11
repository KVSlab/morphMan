import math
from os import makedirs

import vtk
from scipy.interpolate import splrep, splev
from scipy.signal import resample
from vtk.util import numpy_support

from morphman.common.common import *
from morphman.common.vessel_reconstruction_tools import create_parent_artery_patches
from morphman.common.vmtk_wrapper import *
from morphman.common.vmtkpointselector import vmtkPickPointSeedSelector
from morphman.common.vtk_wrapper import *


######################
# SURFACE OPERATIONS #
######################

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
    if not path.isfile(base_path + "_info.txt"):
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
        area.append(vtk_compute_surface_area(delaunay.GetOutput()))
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
        dir_path (str): Location of into.txt file

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
    if flowext or not path.isfile(base_path + "_info.txt"):
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

    cell_normals = vtk_compute_polydata_normals(surface)

    gradients = vtk_compute_normal_gradients(cell_normals)

    # Compute the magnitude of the gradient
    gradients_array = get_vtk_array("Gradients", gradients, 9)
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
        area.append(vtk_compute_surface_area(regions[-1]))

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
                           test_merge=False, changed=False, old_centerline=None):
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
        locator = vtk_point_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(center))
        distances.append(get_distance(region_point, center))

    # Remove the region with the closest distance
    regions.pop(distances.index(min(distances)))

    # Add the other regions back to the surface
    surface = vtk_append_polydata(regions + [surface])
    surface = vtk_clean_polydata(surface)
    surface = vtk_triangulate_surface(surface)

    return surface


def prepare_voronoi_diagram(capped_surface, centerlines, base_path, smooth, smooth_factor, no_smooth, no_smooth_point,
                            voronoi, pole_ids):
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

    Returns:
        voronoi (vtkPolyData): Voronoi diagram of surface.
    """
    # Check if a region should not be smoothed
    if smooth and no_smooth:
        no_smooth_path = base_path + "_centerline_no_smooth.vtp"

        if not path.exists(no_smooth_path):
            # Get inlet and outlets
            tol = get_centerline_tolerance(centerlines)
            inlet = extract_single_line(centerlines, 0)
            inlet = inlet.GetPoint(0)  # inlet.GetNumberOfPoints() - 1)
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
                locator = vtk_point_locator(capped_surface)
                for i in range(len(no_smooth_point) // 3):
                    tmp_id = locator.FindClosestPoint(no_smooth_point[3 * i:3 * (i + 1)])
                    parameters["no_smooth_point_{}".format(i)] = capped_surface.GetPoint(tmp_id)
                    outlets += capped_surface.GetPoint(tmp_id)

            # Store parameters
            write_parameters(parameters)

            # Create the centerline
            no_smooth_centerlines, _, _ = compute_centerlines(inlet, outlets, None, capped_surface, resampling=0.1,
                                                              smooth=False, voronoi=voronoi, pole_ids=pole_ids)

            # Extract the centerline region which diverges from the existing centerlines
            no_smooth_segments = []
            for i in range(no_smooth_centerlines.GetNumberOfLines()):
                tmp_line = extract_single_line(no_smooth_centerlines, i)
                div_ids = []
                for j in range(centerlines.GetNumberOfLines()):
                    div_ids.append(get_diverging_point_id(tmp_line, extract_single_line(centerlines, j), tol))
                div_id = max(div_ids)
                no_smooth_segments.append(extract_single_line(tmp_line, 0, start_id=div_id))

            no_smooth_cl = vtk_append_polydata(no_smooth_segments)
            write_polydata(no_smooth_cl, no_smooth_path)
        else:
            no_smooth_cl = read_polydata(no_smooth_path)

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
        if connected_surface.GetNumberOfPoints() != surface.GetNumberOfPoints():
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


######################
# VORONOI OPERATIONS #
######################


def remove_distant_voronoi_points(voronoi, centerline):
    """Take a voronoi diagram and a centerline remove points that are far away.

    Args:
        voronoi (vtkPolyData): Voronoi data.
        centerline (vtkPolyData): centerline.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram without the extreme points
    """
    n = voronoi.GetNumberOfPoints()
    new_voronoi = vtk.vtkPolyData()
    cell_array = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    radius = np.zeros(n)

    locator = vtk_point_locator(centerline)
    radius_array_data = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1
    limit = radius_array_data(0)
    limit = limit * 10

    count = 0
    for i in range(n):
        point = voronoi.GetPoint(i)
        cl_point_id = locator.FindClosestPoint(point)
        cl_point = centerline.GetPoint(cl_point_id)
        dist = get_distance(point, cl_point)
        if dist / 3 > radius_array_data(i) or radius_array_data(i) > limit:
            count += 1
            continue

        points.InsertNextPoint(point)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i - count)
        value = radius_array_data(i)
        radius[i - count] = value

    print("Removed %s points from the voronoi diagram" % count)

    radius_array = get_vtk_array(radiusArrayName, 1, n - count)
    for i in range(n - count):
        radius_array.SetTuple(i, [float(radius[i])])

    new_voronoi.SetPoints(points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)

    return new_voronoi


def smooth_voronoi_diagram(voronoi, centerlines, smoothing_factor, no_smooth_cl=None):
    """
    Smooth voronoi diagram based on a given
    smoothingfactor. Each voronoi point
    that has a radius less then MISR*(1-smoothingFactor)
    at the closest centerline point is removed.

    Args:
        voronoi (vtkPolyData): Voronoi diagram to be smoothed.
        centerlines (vtkPolyData): Centerline data.
        smoothing_factor (float): Smoothing factor.
        no_smooth_cl (vktPolyData): Unsmoothed centerline.

    Returns: smoothedDiagram (vtkPolyData): Smoothed voronoi diagram.
    """
    number_of_points = voronoi.GetNumberOfPoints()
    thresholds = get_point_data_array(radiusArrayName, centerlines) * (1 - smoothing_factor)

    # Do not smooth inlet and outlets, set threshold to -1
    start = 0
    end = 0
    for i in range(centerlines.GetNumberOfLines()):
        line = extract_single_line(centerlines, i)
        length = get_curvilinear_coordinate(line)
        end_ = line.GetNumberOfPoints() - 1
        end += end_

        # Point buffer start
        end_id = end_ - np.argmin(np.abs(-(length - length.max()) - thresholds[end]))
        start_id = np.argmin(np.abs(length - thresholds[start]))

        thresholds[start:start + start_id] = -1
        thresholds[end - end_id:end] = -1
        start += end_ + 1
        end += 1

    locator = vtk_point_locator(centerlines)
    if no_smooth_cl is not None:
        no_locator = vtk_point_locator(no_smooth_cl)

    smoothed_diagram = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array_numpy = np.zeros(number_of_points)

    count = 0
    for i in range(number_of_points):
        point = voronoi.GetPoint(i)
        radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        id_ = locator.FindClosestPoint(point)
        cl_point = centerlines.GetPoint(id_)

        if get_distance(point, cl_point) > 2 * thresholds[id_] / (1 - smoothing_factor):
            points.InsertNextPoint(point)
            cell_array.InsertNextCell(1)
            cell_array.InsertCellPoint(count)
            radius_array_numpy[count] = radius
            count += 1

        elif no_smooth_cl is not None:
            dist1 = get_distance(point, centerlines.GetPoint(id_))
            id_1 = no_locator.FindClosestPoint(point)
            dist2 = get_distance(point, no_smooth_cl.GetPoint(id_1))

            if dist2 < dist1:
                points.InsertNextPoint(point)
                cell_array.InsertNextCell(1)
                cell_array.InsertCellPoint(count)
                radius_array_numpy[count] = radius
                count += 1
            else:
                if radius >= thresholds[id_]:
                    points.InsertNextPoint(point)
                    cell_array.InsertNextCell(1)
                    cell_array.InsertCellPoint(count)
                    radius_array_numpy[count] = radius
                    count += 1
        else:
            if radius >= thresholds[id_]:
                points.InsertNextPoint(point)
                cell_array.InsertNextCell(1)
                cell_array.InsertCellPoint(count)
                radius_array_numpy[count] = radius
                count += 1

        radius_array = get_vtk_array(radiusArrayName, 1, count)

    for i in range(count):
        radius_array.SetTuple1(i, radius_array_numpy[i])

    smoothed_diagram.SetPoints(points)
    smoothed_diagram.SetVerts(cell_array)
    smoothed_diagram.GetPointData().AddArray(radius_array)

    return smoothed_diagram


def create_new_surface(complete_voronoi_diagram, poly_ball_size=[120, 120, 120]):
    """
    Envelops an input voronoi diagram
    into a new surface model at a
    given resolution determined by
    the poly_ball_size.

    Args:
        complete_voronoi_diagram (vtkPolyData): Voronoi diagram
        poly_ball_size (list): List of dimensional resolution of output model

    Returns:
        envelope (vtkPolyData): Enveloped surface model.
    """
    modeller = vmtk_polyball_modeller(complete_voronoi_diagram, poly_ball_size)

    # Write the new surface
    marching_cube = vtk_marching_cube(modeller)
    envelope = marching_cube.GetOutput()

    return envelope


def get_split_voronoi_diagram(voronoi, centerlines):
    """Given two centerlines, and a Voronoi diagram, return two Voronoi diagrams based on
    the distance of the two centerlines.

    Args:
        voronoi (vtkPolyData): Input Voronoi diagram
        centerlines (list): A list of centerlines (vtkPolyData). An entery could
                            alternativly be None as well, the corresponding voronoi
                            diagram would then be None as well.

    Returns
        voronoi2 (list): A list of Voronoi diagrams closest to each centerline.
    """
    n = len(centerlines)
    centerline1 = [centerlines[i] for i in range(n) if centerlines[i] is not None]
    voronoi1 = [vtk.vtkPolyData() for i in range(n) if centerlines[i] is not None]
    points1 = [vtk.vtkPoints() for i in range(n) if centerlines[i] is not None]
    cell_array1 = [vtk.vtkCellArray() for i in range(n) if centerlines[i] is not None]
    radius1 = [np.zeros(voronoi.GetNumberOfPoints()) for i in range(n) if centerlines[i] is not None]
    loc1 = [vtk_point_locator(centerlines[i]) for i in range(n) if centerlines[i] is not None]
    count1 = [0 for i in range(n) if centerlines[i] is not None]

    n1 = len(centerline1)

    get_radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1

    for i in range(voronoi.GetNumberOfPoints()):
        dists = []
        point = voronoi.GetPoint(i)
        radius = get_radius(i)
        for i in range(n1):
            dists.append(get_distance(centerline1[i].GetPoint(loc1[i].FindClosestPoint(point)), point))

        index = dists.index(min(dists))

        points1[index].InsertNextPoint(point)
        radius1[index][count1[index]] = radius
        cell_array1[index].InsertNextCell(1)
        cell_array1[index].InsertCellPoint(count1[index])
        count1[index] += 1

    for i in range(n1):
        voronoi1[i].SetPoints(points1[i])
        voronoi1[i].SetVerts(cell_array1[i])
        tmp_radius1 = create_vtk_array(radius1[i][radius1[i] > 0], radiusArrayName)
        voronoi1[i].GetPointData().AddArray(tmp_radius1)

    if n1 != n:
        voronoi2 = []
        for i in range(n):
            if centerlines[i] is None:
                voronoi2.append(None)
            else:
                voronoi2.append(voronoi1[i])
    else:
        voronoi2 = voronoi1

    return voronoi2


#########################
# CENTERLINE OPERATIONS #
#########################

def extract_ica_centerline(base_path, resampling_step, relevant_outlets=None):
    """
    Extract a centerline from the inlet to the first branch.

    Args:
        base_path (str): Path to the case folder.
        resampling_step (float): Resampling step length of the extracted centerline.
        relevant_outlets (ndarray): Array containing points corresponding to outlets

    Returns:
        centerline (vtkPolyData): Extracted centerline.
    """
    # TODO: Extract ICA centerline by comparing daughter branch cross-section areas
    centerlines_path = base_path + "_centerline.vtp"
    input_filepath = base_path + ".vtp"
    ica_centerline_path = base_path + "_ica.vtp"
    centerline_relevant_outlets_path = base_path + "_centerline_relevant_outlets_landmark.vtp"

    if path.exists(ica_centerline_path):
        return read_polydata(ica_centerline_path)

    # Prepare surface and identify in/outlets
    surface, capped_surface = prepare_surface(base_path, input_filepath)
    inlet, outlets = get_inlet_and_outlet_centers(surface, base_path)

    if relevant_outlets is not None:
        outlet1, outlet2 = relevant_outlets[:3], relevant_outlets[3:]
        surface_locator = vtk_point_locator(capped_surface)
        id1 = surface_locator.FindClosestPoint(outlet1)
        id2 = surface_locator.FindClosestPoint(outlet2)
        outlet1 = capped_surface.GetPoint(id1)
        outlet2 = capped_surface.GetPoint(id2)
    else:
        outlet1, outlet2 = get_relevant_outlets(capped_surface, base_path)
    outlets, outlet1, outlet2 = get_sorted_outlets(outlets, outlet1, outlet2, base_path)

    # Compute / import centerlines
    centerlines, _, _ = compute_centerlines(inlet, outlets, centerlines_path,
                                            capped_surface,
                                            resampling=resampling_step,
                                            smooth=False, base_path=base_path)

    # Get relevant centerlines
    centerline_relevant_outlets = compute_centerlines(inlet, outlet1 + outlet2,
                                                      centerline_relevant_outlets_path,
                                                      capped_surface,
                                                      resampling=resampling_step)

    # Extract ICA centerline
    tmp_line_1 = extract_single_line(centerline_relevant_outlets[0], 0)
    tmp_line_2 = extract_single_line(centerline_relevant_outlets[0], 1)
    tolerance = get_centerline_tolerance(tmp_line_1)
    line = extract_single_line(tmp_line_1, 0, start_id=0,
                               end_id=get_diverging_point_id(tmp_line_1, tmp_line_2, tolerance))
    write_polydata(line, ica_centerline_path)

    return line


def get_bifurcating_and_diverging_point_data(centerline, centerline_bif, tol):
    """
    Locate bifurcating point and diverging points
    in a bifurcation.
    End points are set based on the MISR at
    the selected points.

    Args:
        centerline (vtkPolyData): Centerline from inlet to relevant outlets.
        centerline_bif (vtkPolyData): Centerline through bifurcation.
        tol (float): Tolerance parameter.

    Returns:
        data (dict): Contains info about diverging point locations.
    """
    # Sort centerline to start at inlet
    cl1 = extract_single_line(centerline, 0)
    cl2 = extract_single_line(centerline, 1)

    # Declear dictionary to hold results
    data = {"bif": {}, 0: {}, 1: {}}

    # Find lower clipping point
    n_points = min(cl1.GetNumberOfPoints(), cl2.GetNumberOfPoints())
    for i in range(0, n_points):
        point_0 = cl1.GetPoint(i)
        point_1 = cl2.GetPoint(i)
        distance_between_points = get_distance(point_0, point_1)
        if distance_between_points > tol:
            center = cl1.GetPoint(i)
            r = cl1.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
            break

    end, r_end, id_end = move_past_sphere(cl1, center, r, i, step=-1)
    data["bif"]["end_point"] = end
    data["bif"]["div_point"] = center

    # Find the diverging points for the bifurcation
    # continue further downstream in each direction and stop when
    # a point is closer than tol, then move point MISR * X
    locator = vtk_point_locator(centerline_bif)

    for counter, cl in enumerate([cl1, cl2]):
        for i in range(i, cl.GetNumberOfPoints(), 1):
            tmp_point = cl.GetPoint(i)
            closest_point_id = locator.FindClosestPoint(tmp_point)
            closest_point = centerline_bif.GetPoint(closest_point_id)
            distance_between_points = get_distance(tmp_point, closest_point)
            if distance_between_points < tol * 4:
                center = cl.GetPoint(i)
                r = cl1.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
                break

        end, r_end, id_end = move_past_sphere(cl, center, r, i, step=1,
                                              stop=i * 100, scale_factor=1)
        data[counter]["end_point"] = end
        data[counter]["div_point"] = center

    return data


def get_manipulated_centerlines(patch_cl, dx, p1, p2, diverging_id, diverging_centerlines, direction, merge_lines=True):
    """Given a centerline (patch_cl), move the centerline a distance (dx) between two
    points (p1 and p2).

    Args:
        patch_cl (vtkPolyData): Centerlines excluding diverging centerlines.
        dx (ndarray): Direction to move geometry.
        p1 (vtkPolyData): First region point.
        p2: (vtkPolyData): Second region point.
        diverging_id (int): List of index where centerlines diverge from region of interest.
        diverging_centerlines (vtkPolyData): Centerlines which diverge from region of interest.
        direction (str): Manipulation direction parameter.
        merge_lines (bool): Merge centerlines and diverging centerlines.

    Returns:
        centerline (vtkPolyData): Manipulated centerline.
    """
    if diverging_id is not None and merge_lines:
        patch_cl = vtk_append_polydata([patch_cl, diverging_centerlines])

    number_of_points = patch_cl.GetNumberOfPoints()
    number_of_cells = patch_cl.GetNumberOfCells()

    centerline = vtk.vtkPolyData()
    centerline_points = vtk.vtkPoints()
    centerline_cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, number_of_points)

    count = 0
    for i in range(number_of_cells):
        line = extract_single_line(patch_cl, i)
        centerline_cell_array.InsertNextCell(line.GetNumberOfPoints())

        radius_array_data = line.GetPointData().GetArray(radiusArrayName).GetTuple1

        locator = vtk_point_locator(line)
        id1 = locator.FindClosestPoint(p1)
        if diverging_id is not None and i == (number_of_cells - 1):
            # Note: Reuse id2 and idmid from previous loop
            pass
        else:
            id2 = locator.FindClosestPoint(p2)
            idmid = int((id1 + id2) * 0.5)

        for p in range(line.GetNumberOfPoints()):
            point = line.GetPoint(p)
            cl_id = locator.FindClosestPoint(point)

            if direction == "horizont":
                if diverging_id is not None and i == (number_of_cells - 1) and diverging_id < cl_id:
                    dist = -dx * (diverging_id - idmid) ** 0.5 / (diverging_id - idmid) ** 0.5
                else:
                    if cl_id < id1:
                        dist = dx
                    elif id1 <= cl_id < idmid:
                        dist = dx * (idmid ** 2 - cl_id ** 2) / (idmid ** 2 - id1 ** 2)
                    elif idmid <= cl_id < (id2 - 1):
                        dist = -dx * (cl_id - idmid) ** 0.5 / (id2 - idmid) ** 0.5
                    else:
                        dist = -dx

            elif direction == "vertical":
                if diverging_id is not None and i == (number_of_cells - 1) and diverging_id < cl_id:
                    dist = 4 * dx * (diverging_id - id1) * (id2 - diverging_id) / (id2 - id1) ** 2
                else:
                    if id1 <= cl_id <= id2:
                        dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2
                    else:
                        dist = 0

            point = np.asarray(point)
            centerline_points.InsertNextPoint(point + dist)
            radius_array.SetTuple1(count, radius_array_data(p))
            centerline_cell_array.InsertCellPoint(count)
            count += 1

    centerline.SetPoints(centerline_points)
    centerline.SetLines(centerline_cell_array)
    centerline.GetPointData().AddArray(radius_array)

    return centerline


def get_centerline_between_clipping_points(centerline_relevant_outlets, data):
    """Get the centerline between two clipping points.

    Args:
        centerline_relevant_outlets (vtkPolyData): Centerline to the two relevant outlets.
        data (dict): A dictionary data.

    Returns:
        centerline (vtkPolyData): Return clipped centerline.
    """
    line0 = extract_single_line(centerline_relevant_outlets, 0)
    line1 = extract_single_line(centerline_relevant_outlets, 1)
    lines = []
    for i, line in enumerate([line0, line1]):
        loc = vtk_point_locator(line)
        tmp_id_dau = loc.FindClosestPoint(data[i]["end_point"])
        tmp_id_bif = loc.FindClosestPoint(data["bif"]["end_point"])
        lines.append(extract_single_line(line, 0, start_id=tmp_id_bif, end_id=tmp_id_dau))

    centerline = vtk_append_polydata(lines)

    return centerline


def get_diverging_point_id(centerline1, centerline2, tol):
    """
    Find ID of diverging point;
    where two input centerlines diverge
    due to a bifurcation.

    Args:
        centerline1 (vtkPolyData): First centerline.
        centerline2 (vtkPolyData): Second centerline.
        tol (float): Tolerance.

    Returns:
        i (int): ID at diverging point.
    """
    # Find clipping points
    n_points = min(centerline1.GetNumberOfPoints(), centerline2.GetNumberOfPoints())
    get_point1 = centerline1.GetPoints().GetPoint
    get_point2 = centerline2.GetPoints().GetPoint

    for i in range(0, n_points):
        distance_between_points = get_distance(get_point1(i), get_point2(i))
        if distance_between_points > tol:
            break

    return i


def get_curvilinear_coordinate(line):
    """
    Get curvilinear coordinates along
    an input centerline.

    Args:
        line (vtkPolyData): Input centerline

    Returns:
        curv_coor (ndarray): Array of abscissa points.
    """
    curv_coor = np.zeros(line.GetNumberOfPoints())
    for i in range(line.GetNumberOfPoints() - 1):
        pnt1 = np.asarray(line.GetPoints().GetPoint(i))
        pnt2 = np.asarray(line.GetPoints().GetPoint(i + 1))
        curv_coor[i + 1] = np.sqrt(np.sum((pnt1 - pnt2) ** 2)) + curv_coor[i]

    return curv_coor


def get_centerline_tolerance(centerline, n=50):
    """
    Finds tolerance based on
    average length between first N points
    along the input centerline.

    Args:
        centerline (vtkPolyData): Centerline data.
        n (int): Number of points

    Returns:
        tolerance (float): Tolerance value.
    """
    line = extract_single_line(centerline, 0)
    length = get_curvilinear_coordinate(line)
    tolerance = np.mean(length[1:n] - length[:n - 1]) / divergingRatioToSpacingTolerance

    return tolerance


def get_clipped_diverging_centerline(centerline, clip_start_point, clip_end_id):
    """
    Clip the opthamlic artery if present.

    Args:
        centerline (vtkPolyData): Line representing the opthalmic artery centerline.
        clip_start_point (tuple): Point at entrance of opthalmic artery.
        clip_end_id (int): ID of point at end of opthalmic artery.

    Returns:
        patch_eye (vtkPolyData): Voronoi diagram representing opthalmic artery.
    """
    points = [clip_start_point, centerline.GetPoint(clip_end_id)]
    div_points = vtk.vtkPoints()
    for p in points:
        div_points.InsertNextPoint(p)

    patch_cl = create_parent_artery_patches(centerline, div_points, siphon=True)

    return patch_cl


def get_line_to_change(surface, centerline, region_of_interest, method, region_points,
                       stenosis_length):
    """
    Extract and spline part of centerline
    within the geometry where
    area variations will be performed.

    Args:
        surface (vtkPolyData): Surface model.
        centerline (vtkPolyData): Centerline in geometry.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'first_line']
        method (str): Determines which kind of manipulation is performed.
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint.
        stenosis_length (float): Multiplier used to determine the length of the stenosis-affected area.

    Returns:
        line_to_change (vtkPolyData): Part of centerline.
    """
    if region_of_interest == "first_line":
        tol = get_centerline_tolerance(centerline)
        line2 = extract_single_line(centerline, 0)
        number_of_points2 = line2.GetNumberOfPoints()

        # Iterate through lines and find diverging point
        n = centerline.GetNumberOfLines()
        point_ids = []
        for j in range(1, n):
            line1 = extract_single_line(centerline, j)
            number_of_points1 = line1.GetNumberOfPoints()

            n_points = min(number_of_points1, number_of_points2)
            for i in range(n_points):
                point1 = line1.GetPoints().GetPoint(i)
                point2 = line2.GetPoints().GetPoint(i)
                if get_distance(point1, point2) > tol:
                    point_id = i
                    break
            point_ids.append(point_id)

        start_id = 0
        end_id = min(point_ids)
        cl_id = point_ids.index(end_id)

        region_points = list(line2.GetPoint(start_id)) + list(line2.GetPoint(end_id))

    elif region_of_interest in ["commandline", "landmarking", "manual"]:
        # Get points from the user
        if region_of_interest == "manual":
            print("\nPlease select region of interest in the render window.")
            stenosis_point_id = vtk.vtkIdList()
            first = True
            while stenosis_point_id.GetNumberOfIds() not in [1, 2]:
                if not first:
                    print("Please provide only one or two points, try again")

                # Select point on surface
                seed_selector = vmtkPickPointSeedSelector()
                seed_selector.SetSurface(surface)
                if method == "variation" or method == "area":
                    seed_selector.text = "Press space to select the start and endpoint of the" + \
                                         " region of interest, 'u' to undo.\n"
                elif method == "stenosis":
                    seed_selector.text = "Press space to select, the center of a new" + \
                                         " stenosis (one point),\nOR place two points on each side" + \
                                         " of an existing stenosis to remove it, \'u\' to undo."
                elif method == "bend":
                    seed_selector.text = "Press space to select the start and end of the" + \
                                         " bend that you want to manipulate, 'u' to undo.\n"
                seed_selector.Execute()
                stenosis_point_id = seed_selector.GetTargetSeedIds()
                first = False

            region_points = []
            for i in range(stenosis_point_id.GetNumberOfIds()):
                region_points += surface.GetPoint(stenosis_point_id.GetId(i))

            print("-- The chosen region points:", region_points)

        # Get locator
        locator = vtk_point_locator(centerline)

        if len(region_points) == 3:
            # Project point onto centerline
            region_points = centerline.GetPoint(locator.FindClosestPoint(region_points))
            point1 = region_points

            # Get relevant line
            tol = get_centerline_tolerance(centerline)
            cl_id = -1
            dist = 1e10
            while dist > tol / 10:
                cl_id += 1
                line = extract_single_line(centerline, cl_id)
                tmp_loc = vtk_point_locator(line)
                tmp_id = tmp_loc.FindClosestPoint(point1)
                dist = get_distance(point1, line.GetPoint(tmp_id))

            # Get length of stenosis
            misr = get_point_data_array(radiusArrayName, line)
            length = stenosis_length * misr[tmp_loc.FindClosestPoint(point1)]

            # Get ids of start and stop
            centerline_length = get_curvilinear_coordinate(line)
            center = centerline_length[tmp_id]
            region_of_interest_id = (center - length <= centerline_length) * (centerline_length <= center + length)
            start_id = np.argmax(region_of_interest_id)
            end_id = region_of_interest_id.shape[0] - 1 - np.argmax(region_of_interest_id[::-1])

        else:
            point1 = region_points[:3]
            point2 = region_points[3:]
            point1 = centerline.GetPoint(locator.FindClosestPoint(point1))
            point2 = centerline.GetPoint(locator.FindClosestPoint(point2))
            region_points[:3] = point1
            region_points[3:] = point2

            distance1 = []
            distance2 = []
            ids1 = []
            ids2 = []
            for i in range(centerline.GetNumberOfLines()):
                line = extract_single_line(centerline, i)
                tmp_loc = vtk_point_locator(line)
                ids1.append(tmp_loc.FindClosestPoint(point1))
                ids2.append(tmp_loc.FindClosestPoint(point2))
                cl_point1 = line.GetPoint(ids1[-1])
                cl_point2 = line.GetPoint(ids2[-1])
                distance1.append(get_distance(point1, cl_point1))
                distance2.append(get_distance(point2, cl_point2))

            tol = get_centerline_tolerance(centerline) / 10
            total_distance = (np.array(distance1) < tol) * (np.array(distance2) < tol)
            cl_id = np.argmax(total_distance)

            if total_distance[cl_id] == 0:
                raise RuntimeError("The two points provided have to be on the same " +
                                   " line (from inlet to outlet), and not at two different" +
                                   " outlets")
            start_id = min(ids1[cl_id], ids2[cl_id])
            end_id = max(ids1[cl_id], ids2[cl_id])

    # Extract and spline a single line
    line_to_change = extract_single_line(centerline, cl_id, start_id=start_id, end_id=end_id)

    remaining_centerlines = []
    diverging_centerlines = []
    start_point = line_to_change.GetPoint(0)
    end_point = line_to_change.GetPoint(line_to_change.GetNumberOfPoints() - 1)
    tol = get_centerline_tolerance(centerline) * 4
    for i in range(centerline.GetNumberOfLines()):
        line = extract_single_line(centerline, i)
        locator = vtk_point_locator(line)
        id1 = locator.FindClosestPoint(start_point)
        id2 = locator.FindClosestPoint(end_point)
        p1_tmp = line.GetPoint(id1)
        p2_tmp = line.GetPoint(id2)
        close_start = get_distance(start_point, p1_tmp) < tol
        close_end = get_distance(end_point, p2_tmp) < tol

        # Check if the centerline is going through both or None of the start and end of
        # the region of interest
        if close_start == close_end:
            if close_start:
                if start_id != 0:
                    tmp = extract_single_line(centerline, i, start_id=0, end_id=start_id - 1)
                    remaining_centerlines.append(tmp)
                tmp = extract_single_line(centerline, i, start_id=end_id + 1,
                                          end_id=line.GetNumberOfPoints() - 1)
                remaining_centerlines.append(tmp)
            else:
                remaining_centerlines.append(line)
        else:
            diverging_centerlines.append(line)

    # Find diverging ids
    diverging_ids = []
    main_line = extract_single_line(centerline, 0)
    id_start = 0
    for line in diverging_centerlines:
        id_end = min([line.GetNumberOfPoints(), main_line.GetNumberOfPoints()])
        for i in np.arange(id_start, id_end):
            p_div = np.asarray(line.GetPoint(i))
            p_cl = np.asarray(main_line.GetPoint(i))
            if get_distance(p_div, p_cl) > tol * 10:
                diverging_ids.append(i)
                break

    # Spline the single line
    line_to_change = compute_splined_centerline(line_to_change, nknots=25, isline=True)

    if len(diverging_centerlines) == 0:
        diverging_centerlines = None
    remaining_centerlines = vtk_append_polydata(remaining_centerlines)

    return line_to_change, remaining_centerlines, diverging_centerlines, region_points, diverging_ids


def get_region_of_interest_and_diverging_centerlines(centerlines_complete, region_points):
    """Extract the centerline between the region points.

    Args:
        centerlines_complete (vktPolyData): Complete set of centerlines in geometry.
        region_points (ndarray): Two points determining the region of interest.

    Returns:
        centerlines (vtkPolyData): Centerlines excluding divering lines.
        diverging_centerlines (vtkPolyData): Centerlines diverging from the region of interest.
        region_points (ndarray): Sorted region points.
        region_points_vtk (vtkPoints): Sorted region points as vtkData.
        diverging_ids (list): List of indices where a diverging centerline starts.
    """
    centerlines = []
    diverging_centerlines = []
    p1 = region_points[0]
    p2 = region_points[1]

    # Search for divering centerlines
    tol = get_centerline_tolerance(centerlines_complete) * 4
    for i in range(centerlines_complete.GetNumberOfLines()):
        line = extract_single_line(centerlines_complete, i)
        locator = vtk_point_locator(line)
        id1 = locator.FindClosestPoint(p1)
        id2 = locator.FindClosestPoint(p2)
        p1_tmp = line.GetPoint(id1)
        p2_tmp = line.GetPoint(id2)
        if get_distance(p1, p1_tmp) < tol and get_distance(p2, p2_tmp) < tol:
            centerlines.append(line)
        else:
            diverging_centerlines.append(line)

    # Sort and set clipping points to vtk object
    centerline = centerlines[0]
    locator = vtk_point_locator(centerline)
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    if id1 > id2:
        region_points = region_points[::-1]
        id1, id2 = id2, id1

    region_points_vtk = vtk.vtkPoints()
    for point in np.asarray(region_points):
        region_points_vtk.InsertNextPoint(point)

    # Find diverging point(s)
    diverging_ids = []
    for line in diverging_centerlines:
        id_end = min([line.GetNumberOfPoints(), centerline.GetNumberOfPoints()])
        for i in np.arange(id1, id_end):
            p_div = np.asarray(line.GetPoint(i))
            p_cl = np.asarray(centerline.GetPoint(i))
            if get_distance(p_div, p_cl) > tol:
                diverging_ids.append(i)
                break

    centerlines = vtk_append_polydata(centerlines)
    diverging_centerlines = vtk_append_polydata(diverging_centerlines) if len(diverging_centerlines) > 0 else None

    return centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids


def compute_discrete_derivatives(line, neigh=10):
    """Compute the curvature and torsion of a line using 'neigh' number of neighboring points.

    Args:
        line (vtkPolyData): Line to compute geometry from.
        neigh (int): Number of naboring points.

    Returns:
        line (vtkPolyData): Output line with geometrical parameters.
        curv (vtkPolyData): Output line with geometrical parameters.
    """
    n_points = line.GetNumberOfPoints()

    # Compute cumulative chord length
    t = np.zeros(n_points)
    p = []
    for i in range(n_points):
        p.append(np.array(list(line.GetPoint(i))))
        p[i] = np.array(p[i])

    norms = [la.norm(p[j] - p[j - 1]) for j in range(1, n_points)]
    s = sum(norms)
    for i in range(1, n_points):
        s1 = sum(norms[:i + 1])
        t[i] = s1 / s

    # Radius of sliding neighbourhood
    m = neigh

    dxdt = np.zeros(n_points)
    dydt = np.zeros(n_points)
    dzdt = np.zeros(n_points)

    x = np.zeros(n_points)
    y = np.zeros(n_points)
    z = np.zeros(n_points)

    for i in range(n_points):
        x[i] = p[i][0]
        y[i] = p[i][1]
        z[i] = p[i][2]

    for i in range(0, m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(0, 2 * m + 1)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(0, 2 * m + 1)]) / t_sum

    for i in range(m, n_points - m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(i - m, i + m + 1)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(i - m, i + m + 1)]) / t_sum

    for i in range(n_points - m, n_points):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(n_points - 2 * m, n_points)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum

    dgammadt = []
    dgammadt_norm = np.zeros(n_points)
    for i in range(n_points):
        dgammadt.append(np.array([dxdt[i], dydt[i], dzdt[i]]))
        dgammadt_norm[i] = la.norm(dgammadt[i])

    tg = []
    for i in range(n_points):
        tg.append(dgammadt[i] / dgammadt_norm[i])

    t1 = np.zeros(n_points)
    t2 = np.zeros(n_points)
    t3 = np.zeros(n_points)

    for i in range(n_points):
        t1[i] = tg[i][0]
        t2[i] = tg[i][1]
        t3[i] = tg[i][2]

    dt1dt = np.zeros(n_points)
    dt2dt = np.zeros(n_points)
    dt3dt = np.zeros(n_points)

    for i in range(0, m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(0, 2 * m + 1)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(0, 2 * m + 1)]) / t_sum

    for i in range(m, n_points - m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(i - m, i + m + 1)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(i - m, i + m + 1)]) / t_sum

    for i in range(n_points - m, n_points):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(n_points - 2 * m, n_points)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(n_points - 2 * m - 1, n_points)]) / t_sum

    dtgdt = []
    dtgdt_norm = np.zeros(n_points)
    for i in range(n_points):
        dtgdt.append(np.array([dt1dt[i], dt2dt[i], dt3dt[i]]))
        dtgdt_norm[i] = la.norm(dtgdt[i])

    curv = np.zeros(n_points)
    for i in range(n_points):
        curv[i] = dtgdt_norm[i] / dgammadt_norm[i]

    curv = resample(curv, n_points)

    return line, curv


def get_k1k2_basis(curvature, line):
    """
    Create a k1-k2 basis used to determine
    the location of each bend of the carotid
    siphon.

    Args:
        curvature (floats): Curvature array.
        line (vtk): Centerline points.

    Returns:
        line (vtk): Centerline points including k1-k2 basis.
    """

    # The ParallelTransportNormals and the FrenetTangent is not orthonormal
    # (but close) from vmtk. Using the GramSchmidt proses gives E2 and fixes
    # the non-orthogonality
    E1 = get_point_data_array("ParallelTransportNormals", line, k=3)
    T = get_point_data_array("FrenetTangent", line, k=3)
    E2 = np.zeros((E1.shape[0], 3))

    for i in range(E1.shape[0]):
        V = np.eye(3)
        V[:, 0] = T[i, :]
        V[:, 1] = E1[i, :]
        V = gram_schmidt(V)

        E1[i, :] = V[:, 1]
        E2[i, :] = V[:, 2]

    # Compute k_1, k_2 furfilling T' = curv(s)*N(s) = k_1(s)*E_1(s) + k_2(s)*E_2(s).
    # This is simply a change of basis for the curvature vector N. The room
    # of k_1 and k_2 can be used to express both the curvature and the
    # torsion.
    N = get_point_data_array("FrenetNormal", line, k=3)

    k2 = (curvature.T * (E1[:, 1] * N[:, 0] - N[:, 1] * E1[:, 0]) /
          (E2[:, 1] * E1[:, 0] - E2[:, 0] * E1[:, 1]))[0]
    k1 = (-(curvature.T * N[:, 0] + k2 * E2[:, 0]) / E1[:, 0])[0]

    for k in [(k1, "k1"), (k2, "k2")]:
        k_array = create_vtk_array(k[0], k[1])
        line.GetPointData().AddArray(k_array)

    return line


def compute_splined_centerline(line, get_curv=False, isline=False, nknots=50, get_stats=True, get_misr=True):
    """
    Given the knots and coefficients of a B-spline representation,
    evaluate the value of the smoothing polynomial and its derivatives.
    This is a wrapper around the FORTRAN routines splev and splder of FITPACK.

    Args:
        line (vtkPolyData): Centerline points.
        get_curv (bool): Computes curvature profile if True.
        isline (bool): Determines if centerline object is a line or points.
        nknots (int): Number of knots.
        get_stats (bool): Determines if curve attribuites are computed or not.
        get_misr (bool): Determines if MISR values are computed or not.

    Returns:
        line (vtkPolyData): Splined centerline data.
    Returns:
        curv (ndarray): Curvature profile.
    """

    if not isline:
        # Create edges between new_centerline points
        pts = vtk.vtkPoints()
        for i in range((line.GetNumberOfCells())):
            pts.InsertNextPoint(line.GetPoint(i))

        lines = vtk.vtkCellArray()
        for i in range((line.GetNumberOfCells()) - 1):
            newline = vtk.vtkLine()
            newline.GetPointIds().SetId(0, i)
            newline.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(newline)

        line = vtk.vtkPolyData()
        line.SetPoints(pts)
        line.SetLines(lines)

    # Collect data from centerline
    data = np.zeros((line.GetNumberOfPoints(), 3))
    if get_misr:
        misr = get_point_data_array(radiusArrayName, line)

    curv_coor = get_curvilinear_coordinate(line)
    for i in range(data.shape[0]):
        data[i, :] = line.GetPoint(i)

    t = np.linspace(curv_coor[0], curv_coor[-1], nknots + 2)[1:-1]
    fx = splrep(curv_coor, data[:, 0], k=4, t=t)
    fy = splrep(curv_coor, data[:, 1], k=4, t=t)
    fz = splrep(curv_coor, data[:, 2], k=4, t=t)

    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    if get_misr:
        data = np.zeros((len(curv_coor), 4))
        data[:, 3] = misr[:, 0]
        header = ["X", "Y", "Z", radiusArrayName]
    else:
        data = np.zeros((len(curv_coor), 3))
        header = ["X", "Y", "Z"]

    data[:, 0] = fx_
    data[:, 1] = fy_
    data[:, 2] = fz_

    line = convert_numpy_data_to_polydata(data, header)

    # Let vmtk compute curve attributes
    if get_stats:
        line = vmtk_compute_geometric_features(line, smooth=False)

    if get_curv:
        # Analytical curvature
        dlsfx = splev(curv_coor, fx, der=1)
        dlsfy = splev(curv_coor, fy, der=1)
        dlsfz = splev(curv_coor, fz, der=1)

        ddlsfx = splev(curv_coor, fx, der=2)
        ddlsfy = splev(curv_coor, fy, der=2)
        ddlsfz = splev(curv_coor, fz, der=2)

        c1xc2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
        c1xc2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
        c1xc2_3 = ddlsfy * dlsfx - ddlsfx * dlsfy

        curvature = np.sqrt(c1xc2_1 ** 2 + c1xc2_2 ** 2 + c1xc2_3 ** 2) / (dlsfx ** 2 + dlsfy ** 2 + dlsfz ** 2) ** 1.5

        return line, curvature
    else:
        return line

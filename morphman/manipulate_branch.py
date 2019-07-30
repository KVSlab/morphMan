##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from scipy.signal import argrelextrema

# Local import
from morphman.common.argparse_common import *
from morphman.common.surface_operations import *

MAX_BIFURCATION_LENGTH = 25


def manipulate_branch(input_filepath, output_filepath, smooth, smooth_factor, poly_ball_size, no_smooth,
                      no_smooth_point, resampling_step, polar_angle, azimuth_angle, remove_branch,
                      branch_to_manipulate_number, branch_location, translation_method, clamp_branch):
    """
    Primary script for moving or removing a selected branch of any blood vessel.
    Relies on a surface geometry of a 3D blood vessel network.

    Defines in- and out-put files, computes voronoi diagram, and
    moves voronoi diagram according to user selected points.
    Used selects two points defining a diverging branch, and the
    one point defining the new location on the surface.
    Proceedes with either removal of a single branch,
    or traslationg and rotation of a single branch.

    Args:
        clamp_branch (bool): Clamps branch at endpoint if True
        remove_branch (bool: If true, removes selected branch completely.
        input_filepath (str): Path to input surface.
        output_filepath (str): Path to output the manipulated surface.
        smooth (bool): Smooth the Voronoi diagram.
        smooth_factor (float): Smoothing factor used for Voronoi diagram smoothing.
        poly_ball_size (list): Resolution of polyballs used to create surface.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
        resampling_step (float): Resampling length when resampling centerline.
        azimuth_angle (float): Angle to rotate around new surface normal vector. Default is no rotation.
        polar_angle (float): Angle to rotate around the axis spanned by cross product of new normal and frenet normal.
        branch_to_manipulate_number (int): Number of branch to manipulate, ordered from up- to down-stream.
        branch_location (ndarray): Point where branch to manipulate will be moved. Closest point on surface is selected.
        translation_method (str): Method for translating the branch to be manipulated.
    """

    # Input filenames
    base_path = get_path_names(input_filepath)
    centerlines_path = base_path + "_centerline.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)
    inlet, outlets = compute_centers(surface, base_path)

    print("-- Computing centerlines and Voronoi diagram")
    centerlines_complete, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path, capped_surface,
                                                                  resampling=resampling_step, smooth=False,
                                                                  base_path=base_path)

    if smooth:
        print("-- Smoothing Voronoi diagram")
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines_complete, base_path, smooth, smooth_factor,
                                          no_smooth, no_smooth_point, voronoi, pole_ids, resampling_step)

    # Select diverging branch from Brancher and compare with diverging branch found manually
    branches_complete = get_all_branches(centerlines_complete)

    # Select branch to manipulate
    if branch_to_manipulate_number is not None:
        check_branch_number(branch_to_manipulate_number, branches_complete)
        branch_to_manipulate = branches_complete[branch_to_manipulate_number - 1]
    else:
        branch_to_manipulate = pick_branch(capped_surface, branches_complete)

    if not remove_branch:
        # Get new position of branch on model surface
        if translation_method == 'manual':
            new_branch_pos_id, new_branch_pos = pick_new_branch_position(capped_surface)
        elif translation_method == 'commandline':
            new_branch_pos_id, new_branch_pos = get_new_branch_position(branch_location, capped_surface)
        elif translation_method == 'no_translation':
            new_branch_pos_id, new_branch_pos = None, None

    branch_to_manipulate_end_point = get_end_point(branch_to_manipulate)
    centerlines = filter_centerlines(centerlines_complete, branch_to_manipulate_end_point)

    print("-- Clipping Voronoi diagram")
    # Get Voronoi of branch
    voronoi_branch, _ = get_split_voronoi_diagram(voronoi, [branch_to_manipulate, centerlines])
    voronoi_branch, _ = filter_voronoi(voronoi_branch, branch_to_manipulate)

    # Get Voronoi of the remaining geometry
    centerline_for_splitting_voronoi = get_centerline_for_splitting_voronoi(centerlines_complete,
                                                                            branch_to_manipulate.GetPoint(0),
                                                                            base_path,
                                                                            capped_surface,
                                                                            voronoi,
                                                                            pole_ids,
                                                                            resampling_step)
    voronoi_remaining, _ = get_split_voronoi_diagram(voronoi, [centerlines,
                                                               centerline_for_splitting_voronoi])

    if remove_branch:
        print("-- Removing branch")
        detach_branch(voronoi_remaining, centerlines, poly_ball_size, surface, output_filepath, centerlines_complete,
                      base_path)
    else:
        move_and_rotate_branch(polar_angle, azimuth_angle, capped_surface, centerlines, centerlines_complete,
                               branch_to_manipulate, new_branch_pos, new_branch_pos_id, output_filepath,
                               poly_ball_size, surface, voronoi_branch, voronoi_remaining,
                               base_path, translation_method, clamp_branch,
                               centerline_for_splitting_voronoi)


def get_centerline_for_splitting_voronoi(centerlines, starting_point, base_path,
                                         capped_surface, voronoi, pole_ids,
                                         resampling_step):
    """
    Get a version of the centerline that extends closer to the 'main' centerline.
    The logic to create a second centerline is to reduce the effect of 'bumb' left after
    moving a branch.

    Args:
        centerlines (vtkPolyData): The complete centerline.
        starting_point (list): The starting point of the branch-clipped centerline.
        base_path (str): Path to the working folder and name of the case.

    Returns:
        clipping_centerline (vtkPolyData): The new centerline for splitting the Voronoi diagram
    """

    lines = []
    lines_main = []
    distances = []
    tol = get_centerline_tolerance(centerlines)
    for i in range(centerlines.GetNumberOfLines()):
        line = extract_single_line(centerlines, i)
        locator = get_vtk_point_locator(line)
        tmp_id = locator.FindClosestPoint(starting_point)
        dist = get_distance(line.GetPoint(tmp_id), starting_point)
        if dist < tol:
            lines.append(line)
        else:
            distances.append(dist)
            lines_main.append(line)

    # Main line for comparing
    main_line = lines_main[distances.index(min(distances))]

    # Create a centerline going the opposite way
    inlet = main_line.GetPoint(main_line.GetNumberOfPoints()-1)
    outlet = list(lines[0].GetPoint(lines[0].GetNumberOfPoints()-1) + lines[0].GetPoint(0))
    reversed_cl, _, _ = compute_centerlines(inlet, outlet,
                                            base_path + "_outlet_to_branch.vtp",
                                            capped_surface, voronoi=voronoi,
                                            pole_ids=pole_ids, resampling=resampling_step,
                                            smooth=False)
    # Extract each line
    reversed_cl_main = extract_single_line(reversed_cl, 1)
    reversed_cl = extract_single_line(reversed_cl, 0)

    for i in range(min([main_line.GetNumberOfPoints(), lines[0].GetNumberOfPoints()])):
        p_main = main_line.GetPoint(i)
        p = lines[0].GetPoint(i)
        if get_distance(p_main, p) > tol * 10:
            new_starting_point = p
            break

    for i in range(min([reversed_cl_main.GetNumberOfPoints(), reversed_cl.GetNumberOfPoints()])):
        p_main = reversed_cl_main.GetPoint(i)
        p = reversed_cl.GetPoint(i)
        if get_distance(p_main, p) > tol * 10:
            new_starting_point_reversed = p
            break

    new_lines = []
    for line in lines:
        loc = get_vtk_point_locator(line)
        start_id = loc.FindClosestPoint(new_starting_point)
        new_lines.append(extract_single_line(line, 0, start_id=start_id))

    loc = get_vtk_point_locator(reversed_cl)
    start_id = loc.FindClosestPoint(new_starting_point_reversed)
    new_lines.append(extract_single_line(reversed_cl, 0, start_id=start_id))

    new_centerline = vtk_merge_polydata(new_lines)

    return new_centerline


def detach_branch(voronoi_remaining, centerlines, poly_ball_size, surface, output_filepath,
                  centerlines_complete, base_path):
    """
    Reconstructs the Voronoi diagram, creating a new surface,
    without the selected branch, utimatley removing it.
    Proceeds by computing the corresponding centerline in the new
    geometries, used for postprocessing, geometric analysis and
    meshing.

    Args:
        centerlines (vtkPolyData): Relevant centerlines in new geometry
        centerlines_complete (vtkPolyData): Complete set of centerlines
        surface (vtkPolyData): Surface model
        voronoi_remaining (vtkPolyData): Voronoi diagram of remaining surface model
        base_path (string): Path to save location
        output_filepath (str): Path to output the manipulated surface.
        poly_ball_size (list): Resolution of polyballs used to create surface.
    """
    new_centerlines_path = base_path + "_centerline_removed_branch.vtp"
    new_voronoi_path = base_path + "_voronoi_removed_branch.vtp"

    # Create new voronoi diagram and new centerlines
    write_polydata(voronoi_remaining, new_voronoi_path)
    write_polydata(centerlines, new_centerlines_path)

    new_surface = create_new_surface(voronoi_remaining, poly_ball_size=poly_ball_size)

    print("-- Smoothing, cleaning, and checking surface")
    new_surface = prepare_output_surface(new_surface, surface, centerlines, output_filepath, test_merge=False,
                                         changed=False, old_centerline=centerlines_complete)

    print("\n-- Writing new surface to {}.".format(output_filepath))
    write_polydata(new_surface, output_filepath)


def move_and_rotate_branch(polar_angle, azimuth_angle, capped_surface, centerlines, centerlines_complete,
                           diverging_centerline_branch, new_branch_pos, new_branch_pos_id, output_filepath,
                           poly_ball_size, surface, voronoi_branch, voronoi_remaining,
                           base_path, method, clamp_branch, centerline_for_splitting_voronoi):
    """
    Moves, and/or performs rotation of the voronoi diagram, then
    reconstructs the Voronoi diagram, creating a new surface.
    Proceeds by computing the corresponding centerline in the new
    geometries, used for postprocessing, geometric analysis and
    meshing.

    Args:
        clamp_branch (bool): Clamps branch at endpoint if True
        method (str): Translation method for selecting new branch position.
        capped_surface (vtkPolyData): Capped surface model
        centerlines (vtkPolyData): Relevant centerlines in new geometry
        centerlines_complete (vtkPolyData): Complete set of centerlines
        diverging_centerline_branch (vtkPolyData): Diverging centerline
        new_branch_pos (vtkPoint): Point where branch is moved
        new_branch_pos_id (int): ID of point where branch is moved
        surface (vtkPolyData): Surface model
        voronoi_branch (vtkPolyData): Voronoi diagram of branch
        voronoi_remaining (vtkPolyData): Voronoi diagram of remaining surface model
        base_path (string): Path to save location
        output_filepath (str): Path to output the manipulated surface.
        poly_ball_size (list): Resolution of polyballs used to create surface.
        azimuth_angle (float): Angle to rotate around new surface normal vector. Default is no rotation.
        polar_angle (float): Angle to rotate around the axis spanned by cross product of new normal and frenet normal.
    """
    if method in ['commandline', 'manual']:
        if polar_angle != 0 or azimuth_angle != 0:
            description = "moved_and_rotated"
        else:
            description = "moved"
    elif method == 'no_translation':
        description = 'rotated'

    new_centerlines_path = base_path + "_centerline_{}.vtp".format(description)
    new_voronoi_path = base_path + "_voronoi_{}.vtp".format(description)

    # Get surface normals and setup of paramters
    old_normal = get_estimated_surface_normal(diverging_centerline_branch)
    manipulated_centerline = diverging_centerline_branch
    manipulated_voronoi = voronoi_branch
    if method == 'no_translation':
        new_normal = old_normal
        origin = np.asarray(diverging_centerline_branch.GetPoint(0))

    # Perform manipulation
    if method in ['commandline', 'manual']:
        print("-- Translating branch")
        new_normal = get_exact_surface_normal(capped_surface, new_branch_pos_id)
        manipulated_voronoi, manipulated_centerline, origin = move_branch(centerlines, manipulated_centerline,
                                                                          new_branch_pos, old_normal, new_normal,
                                                                          manipulated_voronoi, clamp_branch)

    if azimuth_angle != 0:
        print("-- Rotating branch")
        manipulated_voronoi, manipulated_centerline = rotate_branch(azimuth_angle, manipulated_centerline,
                                                                    manipulated_voronoi, origin, new_normal,
                                                                    clamp_branch)
    if polar_angle != 0:
        print("-- Rotating branch")
        rotation_axis = get_rotation_axis(manipulated_centerline, new_normal)
        manipulated_voronoi, manipulated_centerline = rotate_branch(polar_angle, manipulated_centerline,
                                                                    manipulated_voronoi, origin, rotation_axis,
                                                                    clamp_branch)

    # Create new voronoi diagram and new centerlines
    new_voronoi = vtk_merge_polydata([voronoi_remaining, manipulated_voronoi])
    write_polydata(new_voronoi, new_voronoi_path)

    new_centerlines = vtk_merge_polydata([centerlines, manipulated_centerline])
    write_polydata(new_centerlines, new_centerlines_path)

    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    print("-- Smoothing, cleaning, and checking surface")
    new_surface = prepare_output_surface(new_surface, surface, new_centerlines, output_filepath, test_merge=False,
                                         changed=True, old_centerline=centerlines_complete)

    print("\n-- Writing new surface to {}.".format(output_filepath))
    write_polydata(new_surface, output_filepath)


def get_rotation_axis(centerline, normal_vector):
    """
    Compute axis of rotation for Azimuthal rotation

    Args:
        centerline (vtkPolyData): Centerline of branch to manipulate
        normal_vector (vtkPolyData): Surface normal vector at rotation origin

    Returns:
        ndarray: Axis of rotation vector
    """
    centerline = vmtk_compute_geometric_features(centerline, True)
    frenet_normal = centerline.GetPointData().GetArray("FrenetNormal").GetTuple3(0)
    rotation_axis = np.cross(normal_vector, frenet_normal)

    return rotation_axis


def move_branch(centerlines, diverging_centerline_branch, new_branch_pos,
                old_normal, new_normal, voronoi_branch, clamp_branch):
    """
    Translates and rotates the voronoi diagram and centerline from one point
    on the surface model to a selected point, defined by
    new branch position.

    Args:
        clamp_branch (bool): Clamps branch at endpoint if True
        new_normal (ndarray): New surface normal
        old_normal (ndarray): Old surface normal
        centerlines (vtkPolyData): Relevant centerlines in new geometry
        diverging_centerline_branch (vtkPolyData): Diverging centerline
        new_branch_pos (vtkPoint): Point where branch is moved
        voronoi_branch (vtkPolyData): Voronoi diagram of branch

    Returns:
        vtkPolyData: Translated Voronoi diagram
        vtkPolyData: Translated centerline
        ndarray: Origin of new position
    """
    # Define rotation between surface normal vectors
    rotation_axis, surface_normals_angle = get_rotation_axis_and_angle(new_normal, old_normal)
    R = get_rotation_matrix(rotation_axis, surface_normals_angle)

    # Define translation parameters
    dx, origin = get_translation_parameters(centerlines, diverging_centerline_branch, new_branch_pos)

    # Move branch centerline and voronoi diagram
    moved_voronoi = manipulate_voronoi_branch(voronoi_branch, dx, R, origin, diverging_centerline_branch, rotation_axis,
                                              surface_normals_angle, 'translate', clamp_branch)
    moved_centerline = manipulate_centerline_branch(diverging_centerline_branch, origin, R, dx, rotation_axis,
                                                    surface_normals_angle, 'translate', clamp_branch)

    return moved_voronoi, moved_centerline, origin


def rotate_branch(angle, diverging_centerline_branch, voronoi_branch, origin, axis_of_rotation, clamp_branch):
    """
    Perform rotation of the voronoi diagram and the centerline
    around a given axis defined by the input normal vector.

    Args:
        clamp_branch (bool): Clamps branch at endpoint if True
        angle (float): Angle to rotate the branch, in radians
        origin (ndarray): Origin of the centerline
        axis_of_rotation (ndarray): Vector defing axis to rotate around
        diverging_centerline_branch (vtkPolyData): Diverging centerline
        voronoi_branch (vtkPolyData): Voronoi diagram of branch

    Returns:
        vtkPolyData: Rotated Voronoi diagram
        vtkPolyData: Rotated centerline
    """
    # Define rotation around new surface normal
    R = get_rotation_matrix(-axis_of_rotation, angle)

    # Move branch centerline and voronoi diagram
    rotated_voronoi = manipulate_voronoi_branch(voronoi_branch, 0.0, R, origin, diverging_centerline_branch,
                                                axis_of_rotation, angle, 'rotate', clamp_branch)
    rotated_centerline = manipulate_centerline_branch(diverging_centerline_branch, origin, R, 0.0, axis_of_rotation,
                                                      angle, 'rotate', clamp_branch)

    return rotated_voronoi, rotated_centerline


def check_branch_number(branch_to_manipulate_number, branches_complete):
    """
    Check if branch number provided by user is larger than number of centerlines.
    Raises RuntimeError if number exceeds limit.

    Args:
        branch_to_manipulate_number (int): Input number, supplied by user
        branches_complete (list): All centerline branches
    """
    num_lines = len(branches_complete)
    if branch_to_manipulate_number > num_lines:
        raise RuntimeError("\nERROR: Branch number cannot exceed number of valid branches." +
                           " Number of selectable branches for this model is {}.".format(num_lines))


def get_new_branch_position(branch_location, capped_surface):
    """
    Get point on surface closest to branch_location. Returns
    both the point on the surface and it's ID.

    Args:
        branch_location (ndarray): User selected point.
        capped_surface (vtkPolyData): Input surface

    Returns:
        new_branch_pos_id (int): Point closest to branch location ID on surface.
        new_branch_pos (ndarray): Point closest to branch location on surface.
    """
    surface_locator = get_vtk_point_locator(capped_surface)
    new_branch_pos_id = surface_locator.FindClosestPoint(branch_location)
    new_branch_pos = capped_surface.GetPoint(new_branch_pos_id)

    return new_branch_pos_id, new_branch_pos


def pick_new_branch_position(capped_surface):
    """
    Select (by manually chosing) where branch is translated to on the surface.

    Args:
        capped_surface (vtkPolyData): Input surface

    Returns:
        new_branch_pos_id (int): Point closest to branch location ID on surface.
        new_branch_pos (ndarray): Point closest to branch location on surface.
    """
    print("\nPlease select new point to move branch in the render window.")
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(capped_surface)
    seed_selector.text = "Press space to select the point where you want to move the branch to, 'u' to undo.\n"
    seed_selector.Execute()
    new_branch_pos_id = seed_selector.GetTargetSeedIds().GetId(0)
    new_branch_pos = capped_surface.GetPoint(new_branch_pos_id)
    print("-- New branch location is:", np.asarray(new_branch_pos))

    return new_branch_pos_id, new_branch_pos


def pick_branch(capped_surface, centerlines):
    """
    Select (by manually chosing) which branch is translated to on the surface.

    Args:
        capped_surface (vtkPolyData): Input surface
        centerlines (vtkPolyData): Set of valid centerline branches throughout model

    Returns:
          branch_to_manipulate (vtkPolyData): Branch selected to be manipulated
    """
    print("\nPlease select the branch you with to move in the render window.")
    # Select point on surface
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(capped_surface)
    seed_selector.text = "Press space to select the branch you want to move, 'u' to undo.\n"
    seed_selector.Execute()
    branch_surface_point_id = seed_selector.GetTargetSeedIds().GetId(0)
    surface_point = capped_surface.GetPoint(branch_surface_point_id)

    # Find closest centerlines
    closest_dist = 1e9
    valid_branches_to_manipulate = []
    for line_to_compare in centerlines:
        locator = get_vtk_point_locator(line_to_compare)
        closest_point_id = locator.FindClosestPoint(surface_point)
        closest_point = line_to_compare.GetPoint(closest_point_id)
        dist = get_distance(closest_point, surface_point)
        if dist <= closest_dist:
            closest_dist = dist
            valid_branches_to_manipulate.append(line_to_compare)

    # Select shortest centerline of valid branches
    branch_to_manipulate = valid_branches_to_manipulate[0]
    for branch in valid_branches_to_manipulate[1:]:
        branch_length = len(get_curvilinear_coordinate(branch))
        current_length = len(get_curvilinear_coordinate(branch_to_manipulate))
        if branch_length < current_length:
            branch_to_manipulate = branch

    return branch_to_manipulate


def filter_voronoi(voronoi, diverging_centerline_branch):
    """
    Filter away voronoi points too far away from relevant branch

    Args:
        voronoi (vtkPolyData): Voronoi diagram to be filtered
        diverging_centerline_branch (vtkPolyData): Relevant centerlines of the branch

    Returns:
        vtkPolyData: Voronoi diagram, diverging part
    Returns:
        vtkPolyData: Voronoi diagram, remaining part
    """
    misr = get_point_data_array(radiusArrayName, diverging_centerline_branch)
    locator = get_vtk_point_locator(diverging_centerline_branch)

    n = voronoi.GetNumberOfPoints()
    diverging_voronoi = vtk.vtkPolyData()
    remaining_voronoi = vtk.vtkPolyData()
    div_cell_array = vtk.vtkCellArray()
    rem_cell_array = vtk.vtkCellArray()
    div_points = vtk.vtkPoints()
    div_radius = np.zeros(n)
    rem_points = vtk.vtkPoints()
    rem_radius = np.zeros(n)

    radius_array_data = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1

    div_count = 0
    rem_count = 0
    for i in range(n):
        point = voronoi.GetPoint(i)
        closest_point_id = locator.FindClosestPoint(point)
        closest_point = diverging_centerline_branch.GetPoint(closest_point_id)
        if get_distance(closest_point, point) > max(misr):
            rem_count = set_voronoi_point_data(i, point, radius_array_data, rem_cell_array, rem_count, rem_points,
                                               rem_radius)
        else:
            div_count = set_voronoi_point_data(i, point, radius_array_data, div_cell_array, div_count, div_points,
                                               div_radius)

    set_voronoi_data(div_cell_array, div_count, div_points, div_radius, diverging_voronoi)
    set_voronoi_data(rem_cell_array, rem_count, rem_points, rem_radius, remaining_voronoi)

    return diverging_voronoi, remaining_voronoi


def set_voronoi_point_data(i, point, radius_array_data, cell_array, count, points, radius):
    """
    Set point data to a single Voronoi diagram point

    Args:
        i (int): Counter
        point (vtkPoint): Single point of Voronoi diagram
        radius_array_data (ndarray): MISR Radius array
        cell_array (vtkCellArray): Cell array
        count (int): Specific counter
        points (vtkPoints): Point array
        radius (ndarray):  Radius array

    Returns:
        int: Incremented counter
    """
    points.InsertNextPoint(point)
    cell_array.InsertNextCell(1)
    cell_array.InsertCellPoint(count)
    value = radius_array_data(i)
    radius[count] = value
    count += 1

    return count


def set_voronoi_data(cell_array, count, points, radius, voronoi):
    """
    Apply points and data to voronoi object

    Args:
        cell_array (vtkCellArray): Cell array
        count (int): Specific counter
        points (vtkPoints): Point array
        radius (ndarray):  Radius array
        voronoi (vtkPolyData): Voronoi diagram to be filtered
    """
    radius_array = get_vtk_array(radiusArrayName, 1, count)
    for i in range(count):
        radius_array.SetTuple(i, [float(radius[i])])

    voronoi.SetPoints(points)
    voronoi.SetVerts(cell_array)
    voronoi.GetPointData().AddArray(radius_array)


def get_translation_parameters(centerlines, diverging_centerline_branch, new_branch_pos):
    """
    Get distance to translate branch, and the new position location for branch

    Args:
        centerlines (vtkPolyData): Relevant centerlines
        diverging_centerline_branch (vtkPolyData): Diverging centerline
        new_branch_pos (vtkPoint):  Position where branch is moved.

    Returns:
        dx (float): Distance to translate branch
    Returns:
        origin (ndarray): Adjusted origin / position of new branch position
    """
    locator = get_vtk_point_locator(centerlines)
    new_branch_pos_on_cl = centerlines.GetPoint(locator.FindClosestPoint(new_branch_pos))
    adjusted_branch_pos = (np.asarray(new_branch_pos) - np.asarray(new_branch_pos_on_cl)) * 0.8 + np.asarray(
        new_branch_pos_on_cl)
    dx = np.asarray(adjusted_branch_pos) - np.asarray(diverging_centerline_branch.GetPoint(0))
    origin = np.asarray(adjusted_branch_pos)

    return dx, origin


def get_exact_surface_normal(capped_surface, new_branch_pos_id):
    """
    Compute normal out of surface at given point.

    Args:
        capped_surface (vtkPolyData): Capped surface model
        new_branch_pos_id (int): ID of point where branch is moved, on surface

    Returns:
        new_normal (ndarray): Normal vector out of surface
    """
    capped_surface_with_normals = vmtk_compute_surface_normals(capped_surface)
    new_normal = capped_surface_with_normals.GetPointData().GetNormals().GetTuple(new_branch_pos_id)
    new_normal /= la.norm(new_normal)

    return new_normal


def get_estimated_surface_normal(centerline):
    """
    Estimate the surface normal at initial diverging centerline branch.

    Args:
        centerline (vtkPolyData): Diverging centerline to be moved.

    Returns:
        normal_vector (ndarray): Estimated normal vector at diverging centerline
    """
    centerline = vmtk_compute_geometric_features(centerline, True)
    curvature = get_point_data_array("Curvature", centerline)
    first_local_maxima_id = argrelextrema(curvature, np.greater)[0][0]

    # TODO: Generalize choice of end point factor
    factor = 0.4
    start_point = 0
    end_point = int(first_local_maxima_id * factor)
    normal_vector = np.asarray(centerline.GetPoint(end_point)) - np.asarray(
        centerline.GetPoint(start_point))

    normal_vector /= la.norm(normal_vector)
    return normal_vector


def manipulate_centerline_branch(centerline_branch, origin, R, dx, normal, angle, manipulation, clamp_branch):
    """
    Depending on manipulation method, either translates or
    rotates the selected branch, represented as a centerline.

    Args:
        clamp_branch (bool): Clamps branch at endpoint if True
        manipulation (str): Type of manipulation, either 'rotate' or 'translate'
        centerline_branch (vtkPolyData): Centerline branch to be manipulated
        dx (float): Distance to translate branch
        R (ndarray): Rotation matrix, rotation from old to new surface normal or around new surface normal
        origin (ndarray): Adjusted origin / position of new branch position
        angle (float): Angle to rotate in radians
        normal (ndarray): Normal vector at manipulation location

    Returns:
        centerline (vtkPolyData): Manipulated centerline
    """
    # Locator to find closest point on centerline
    number_of_points = centerline_branch.GetNumberOfPoints()

    centerline = vtk.vtkPolyData()
    centerline_points = vtk.vtkPoints()
    centerline_cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, number_of_points)

    centerline_cell_array.InsertNextCell(number_of_points)
    radius_array_data = centerline_branch.GetPointData().GetArray(radiusArrayName).GetTuple1

    # Transition from base to branch for 30% of branch creating a smoother rotation transition
    base_end = int(number_of_points * 0.3)
    centerline_locator = get_vtk_point_locator(centerline_branch)

    for p in range(centerline_branch.GetNumberOfPoints()):
        point = centerline_branch.GetPoint(p)

        if manipulation == 'translate':
            # Translate and rotate branch upright
            if clamp_branch:
                R, point = get_clamped_branch_translation_factors(angle, centerline_locator, dx, normal,
                                                                  number_of_points, point)
            else:
                point = np.asarray(point) + dx

            point = np.dot(R, point - origin) + origin

        elif manipulation == 'rotate':
            # Rotate branch around axis
            if clamp_branch:
                point = get_clamped_branch_rotation_factors(angle, p, number_of_points, normal, origin, point)
            else:
                if p <= base_end:
                    transition_angle = angle * rotation_profile(p, base_end)
                    transition_rotation_matrix = get_rotation_matrix(-normal, transition_angle)
                    point = np.dot(transition_rotation_matrix, point - origin) + origin
                else:
                    point = np.dot(R, point - origin) + origin

        centerline_points.InsertNextPoint(point)
        radius_array.SetTuple1(p, radius_array_data(p))
        centerline_cell_array.InsertCellPoint(p)

    centerline.SetPoints(centerline_points)
    centerline.SetLines(centerline_cell_array)
    centerline.GetPointData().AddArray(radius_array)

    return centerline


def get_clamped_branch_translation_factors(angle, centerline_locator, dx, normal, number_of_points, point):
    """
    Gradually maniulate branch to clamp branch at endpoint,
    when rotating branch.

    Args:
        dx (float): Distance to translate branch
        angle (float): Angle to rotate in radians
        normal (ndarray): Normal vector at manipulation location
        centerline_locator (vtkPointLocator): Locator for centerline
        number_of_points (int): Number of centerline points
        point (vtkPoint): Voronoi point

    Returns:
        ndarray: Rotation matrix
        ndarray: Translation point
    """
    cl_id = centerline_locator.FindClosestPoint(point)
    clamp_factor = clamp_profile(cl_id, number_of_points)
    new_angle = angle * clamp_factor
    R = get_rotation_matrix(normal, new_angle)
    point = np.asarray(point) + dx * clamp_factor

    return R, point


def clamp_profile(centerline_id, number_of_points):
    """
    Profile used for gradually translating a branch to be clamped.
    Currently using a linear profile, ranging from 0 to 1.

    Args:
        centerline_id (int): ID at current centerline point
        number_of_points (int): Number of centerline points

    Returns:
        float: Clamp factor, ranging from 0 to 1
    """
    return (number_of_points - centerline_id) / float(number_of_points)


def get_rotation_axis_and_angle(new_normal, old_normal):
    """
    Compute axis vector and angle between normal vectors (input)

    Args:
        old_normal (ndarray): Normal vector at initial position
        new_normal (ndarray): Normal vector at new position

    Returns:
        u (ndarray): Normal vector corresponding to rotation axis
    Returns:
        angle (float): Angle between normal vectors
    """
    z = np.asarray(old_normal)
    z_prime = np.asarray(new_normal)

    u = np.cross(z, z_prime)
    u /= la.norm(u)

    angle = get_angle(z, z_prime)

    return u, angle


def manipulate_voronoi_branch(voronoi, dx, R, origin, centerline, normal, angle, manipulation, clamp_branch):
    """
    Depending on manipulation method, either translates or
    rotates the selected branch, represented as a Voronoi diagram.

    Args:
        clamp_branch (bool): Clamps branch at endpoint if True
        manipulation (str): Type of manipulation, either 'rotate' or 'translate'
        voronoi (vtkPolyData): Voronoi diagram of surface
        dx (float): Distance to translate branch
        R (ndarray): Rotation matrix, rotation from old to new surface normal or around new surface normal
        origin (ndarray): Adjusted origin / position of new branch position
        angle (float): Angle to rotate in radians
        normal (ndarray): Normal vector at manipulation location
        centerline (vtkPolyData): Centerline of branch to be manipulated

    Returns:
        new_voronoi (vtkPolyData): Manipulated Voronoi diagram
    """
    # Voronoi diagram
    n = voronoi.GetNumberOfPoints()
    new_voronoi = vtk.vtkPolyData()
    voronoi_points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, n)
    misr_array = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1

    # Centerline locator
    m = centerline.GetNumberOfPoints()
    centerline_locator = get_vtk_point_locator(centerline)

    # Transition from base to branch for 30% of branch creating a smoother rotation transition
    base_end = int(m * 0.3)

    # Iterate through Voronoi diagram and manipulate
    for i in range(n):
        point = voronoi.GetPoint(i)
        misr = misr_array(i)

        if manipulation == 'translate':
            # Translate voronoi points
            if clamp_branch:
                R, point = get_clamped_branch_translation_factors(angle, centerline_locator, dx, normal, m, point)
            else:
                point = np.asarray(point) + dx

            point = np.dot(R, point - origin) + origin

        elif manipulation == 'rotate':
            # Rotate voronoi points
            cl_id = centerline_locator.FindClosestPoint(point)

            if clamp_branch:
                point = get_clamped_branch_rotation_factors(angle, cl_id, m, normal, origin, point)
            else:
                if cl_id <= base_end:
                    transition_angle = angle * rotation_profile(cl_id, base_end)
                    transition_rotation_matrix = get_rotation_matrix(-normal, transition_angle)
                    point = np.dot(transition_rotation_matrix, point - origin) + origin
                else:
                    point = np.dot(R, point - origin) + origin

        # Change radius
        radius_array.SetTuple1(i, misr)
        voronoi_points.InsertNextPoint(point)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i)

    new_voronoi.SetPoints(voronoi_points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)

    return new_voronoi


def rotation_profile(centerline_id, number_of_points):
    """
    Profile used for gradually rotating a branch,
    to avoid artifacts at the base.
    Currently using a root function profile, ranging from 0 to 1.

    Args:
        centerline_id (int): ID at current centerline point
        number_of_points (int): Number of centerline points

    Returns:
        float: Clamp factor, ranging from 0 to 1
    """
    return (centerline_id / number_of_points) ** 0.2


def get_clamped_branch_rotation_factors(angle, cl_id, m, axis_of_rotation, origin, point):
    """
    Gradually maniulate branch to clamp branch at endpoint,
    when rotating branch.

    Args:
        angle (float): Angle to rotate in radians
        axis_of_rotation (ndarray): Axis to rotate around
        point (vtkPoint): Voronoi point
        cl_id (int): Current centerline point ID
        m (int): Number of centerline points
        origin (ndarray): Origin to rotate around

    Returns:
        ndarray: Rotated Voronoi point
    """
    transition_angle = angle * ((cl_id / float(m)) ** 0.2 + clamp_profile(cl_id, m) - 1)
    R = get_rotation_matrix(-axis_of_rotation, transition_angle)
    point = np.dot(R, point - origin) + origin

    return point


def get_all_branches(centerlines):
    """
    Extract and combine all branches of the surface model.
    Removes first part of centerlines after combining,
    excluding the bifurcating part.

    Args:
        centerlines (list, ndarray): List containing sorted centerlines

    Returns:
        list: All possible branches in the geometry
    """

    # Get branches
    branched_centerlines = vmtk_compute_branch_extractor(centerlines)

    # Remove first segment from inlet (cannot remove or move this)
    branched_centerlines = vtk_compute_threshold(branched_centerlines, "TractIds",
                                                 threshold_type="between",
                                                 lower=0.1, upper=1e6, source=1)

    # Remove the bifurcation segments
    branched_centerlines = vtk_compute_threshold(branched_centerlines, "Blanking",
                                                 threshold_type="between", lower=-0.1,
                                                 upper=0.1, source=1)

    tract_ids = get_cell_data_array("TractIds", branched_centerlines, k=1)

    n = centerlines.GetNumberOfLines()
    centerline_lines = [extract_single_line(centerlines, i) for i in range(n)]
    locators = [get_vtk_point_locator(centerline_lines[i]) for i in range(n)]
    tol = get_centerline_tolerance(centerlines)

    # Storing each branch
    branches = []

    for i in np.unique(tract_ids):
        lines = vtk_compute_threshold(branched_centerlines, "TractIds",
                                      threshold_type="between", lower=i-0.5,
                                      upper=i+0.5, source=1)

        # Get unique start points
        start_points = []
        for j in range(lines.GetNumberOfLines()):
            tmp = extract_single_line(lines, j).GetPoint(0)
            if tmp not in start_points:
                start_points.append(tmp)

        # Get all upstream branches of each point
        for start_point in start_points:
            point_ids = [loc.FindClosestPoint(start_point) for loc in locators]
            points = [centerline_lines[i].GetPoint(point_ids[i]) for i in range(n)]
            distances = np.array([get_distance(start_point, points[i]) for i in range(n)])
            branch = [extract_single_line(centerline_lines[i], 0, start_id=point_ids[i]) for i in range(n)
                      if distances[i] > tol]
            branches.append(vtk_merge_polydata(branch))

    return branches


def read_command_line_branch(input_path=None, output_path=None):
    """
    Read arguments from commandline and return all values in a dictionary.
    If input_path and output_path are not None, then do not parse command line, but
    only return default values.

    Args:
        input_path (str): Input file path, positional argument with default None.
        output_path (str): Output file path, positional argument with default None.
    """
    # Description of the script
    description = "Moves a selected part of a tubular geometry, " + \
                  "in two (horizontal, vertical) geometry-specific directions. " + \
                  "Magnitude of movement is defined by the parameters alpha and beta" + \
                  "Primary script used for application in blood vessels."

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    # Add common arguments
    required = not (input_path is not None and output_path is not None)
    add_common_arguments(parser, required=required)

    parser.add_argument("-tm", "--translation-method", type=str, default="manual",
                        choices=["manual", "commandline", "no_translation"],
                        help="Defines the method of translation of the branch to be manipulated." +
                             " The parameter provides three options: 'manual', 'commandline' and 'no_translation'. In" +
                             " 'manual' the user will be provided with a visualization of the input surface, and " +
                             "asked to provide the new position of the branch on the surface model." +
                             " If 'commandline' is provided, then '--branch-location'" +
                             " is expected to be provided. Selecting 'no_translation' will " +
                             "result in no translation; any manipulation performed on the " +
                             "branch will happen at the branch's current position. ")
    parser.add_argument('-bl', "--branch-location", nargs="+", type=float, default=None, metavar="branch_location",
                        help="If this parameter is provided, the branch to be manipulated will be moved to the point "
                             "on the surface closest to this point. Example providing the point (1, 5, -1):" +
                             " --branch-loc 1 5 -1")

    # Arguments for rotation
    parser.add_argument('-aa', '--azimuth-angle', type=float, default=0,
                        help="The manipulated branch is rotated an angle 'aa' around the old or new" +
                             " surface normal vector. 'aa' is assumed to be in degrees," +
                             " and not radians. Default is no rotation.", metavar="surface_normal_axis_angle")
    parser.add_argument('-pa', '--polar-angle', type=float, default=0,
                        help="The manipulated branch is rotated an angle 'pa' around the" +
                             " surface tangent vector, constructed by the cross product of the surface normal vector" +
                             " and the Frenet normal vector. 'pa' is assumed to be in degrees," +
                             " and not radians. Default is no rotation.", metavar="surface_tangent_axis_angle")

    # Argument for selecting branch
    parser.add_argument('-bn', '--branch-number', type=int, default=None,
                        help="The number corresponding the branch to manipulate. " +
                             "The branches are ordered from 1 to N, " +
                             "from upstream to downstream, relative to the inlet. " +
                             "If not selected, the user must manually select the branch "
                             "to manipulate. ", metavar="branch_number")

    # Argument for selecting branch
    parser.add_argument('-rb', '--remove-branch', type=str2bool, default=False,
                        help="If True, will remove selected branch and perform no manipulation")

    # Argument for clamping branch when translating
    parser.add_argument('-cb', '--clamp-branch', type=str2bool, default=False,
                        help="If True, will clamp selected branch to branch endpoint")

    # Parse paths to get default values
    if required:
        args = parser.parse_args()
    else:
        args = parser.parse_args(["-i" + input_path, "-o" + output_path])

    if not 0 <= args.azimuth_angle <= 360:
        raise ArgumentTypeError("The azimuth angle is limited to be within [0, 360] degrees, cannot have value" +
                                " {}".format(args.azimuth_angle))

    if not -180 <= args.polar_angle <= 180:
        raise ArgumentTypeError("The polar angle is limited to be within [-180, 180] degrees, cannot have value" +
                                " {}".format(args.polar_angle))

    # Convert from deg to rad and invert rotation if exceeding 180 degrees
    polar_angle_to_radians = args.polar_angle * math.pi / 180

    azimuth_angle_to_radians = args.azimuth_angle * math.pi / 180
    if azimuth_angle_to_radians > np.pi:
        azimuth_angle_to_radians -= 2 * np.pi

    if args.branch_number is not None:
        if args.branch_number < 1:
            raise ValueError("ERROR: Branch number cannot be 0 or negative. Please select a positive number")

    if args.no_smooth_point is not None and len(args.no_smooth_point):
        if len(args.no_smooth_point) % 3 != 0:
            raise ValueError("ERROR: Please provide the no smooth point(s) as a multiple of 3")

    return dict(input_filepath=args.ifile, smooth=args.smooth, smooth_factor=args.smooth_factor,
                output_filepath=args.ofile, poly_ball_size=args.poly_ball_size, no_smooth=args.no_smooth,
                no_smooth_point=args.no_smooth_point, resampling_step=args.resampling_step,
                polar_angle=polar_angle_to_radians, azimuth_angle=azimuth_angle_to_radians,
                clamp_branch=args.clamp_branch, remove_branch=args.remove_branch,
                branch_to_manipulate_number=args.branch_number, branch_location=args.branch_location,
                translation_method=args.translation_method)


def main_branch():
    manipulate_branch(**read_command_line_branch())


if __name__ == "__main__":
    manipulate_branch(**read_command_line_branch())

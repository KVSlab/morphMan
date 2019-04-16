##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

import functools
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from scipy.signal import argrelextrema

# Local import
from common import *


def manipulate_branch(input_filepath, output_filepath, smooth, smooth_factor, poly_ball_size, no_smooth,
                      no_smooth_point, resampling_step, angle, branch_to_manipulate_number, branch_location):
    """
    Primary script for moving a selected branch of any blood vessel.
    Relies on a surface geometry of a 3D blood vessel network.

    Defines in- and out-put files, computes voronoi diagram, and
    moves voronoi diagram according to user selected points.
    Used selects two points defining a diverging branch, and the
    one point defining the new location on the surface.
    Moves, and performs rotation of the voronoi diagram, then
    reconstructs the Voronoi diagram, creating a new surface.
    Proceeds by computing the corresponding centerline in the new
    geometries, used for postprocessing, geometric analysis and
    meshing.

    Args:
        input_filepath (str): Path to input surface.
        output_filepath (str): Path to output the manipulated surface.
        smooth (bool): Smooth the Voronoi diagram.
        smooth_factor (float): Smoothing factor used for Voronoi diagram smoothing.
        poly_ball_size (list): Resolution of polyballs used to create surface.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
        resampling_step (float): Resampling length when resampling centerline.
        angle (float): Angle to rotate around new surface normal vector. Default i no rotation.
        branch_to_manipulate_number (int): Number of branch to manipulate, ordered from up- to down-stream.
        branch_location (ndarray): Point where branch to manipulate will be moved. Closest point on surface is selected.
    """

    # Input and output filenames
    base_path = get_path_names(input_filepath)
    centerlines_path = base_path + "_centerline.vtp"
    new_centerlines_path = base_path + "_centerline_moved_and_rotated.vtp"
    new_voronoi_path = base_path + "_voronoi_moved_and_rotated.vtp"
    unprepared_output_filepath = base_path + "_unprepared_output.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)
    inlet, outlets = get_centers(surface, base_path)

    print("-- Compute centerlines and Voronoi diagram")
    centerlines_complete, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                                  capped_surface, resampling=resampling_step,
                                                                  smooth=False, base_path=base_path)

    if smooth:
        print("-- Smoothing Voronoi diagram")
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines_complete, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)

    # Select branch to manipulate
    if branch_to_manipulate_number > 0:
        check_branch_number(branch_to_manipulate_number, centerlines_complete)
        longest_centerline = get_sorted_lines(centerlines_complete)[-1]
        branch_to_manipulate = get_branch(branch_to_manipulate_number, centerlines_complete, longest_centerline)
    else:
        shortest_centerlines = get_sorted_lines(centerlines_complete)[:-1]
        branch_to_manipulate = pick_branch(capped_surface, shortest_centerlines)

    # Get new position of branch on model surface
    if branch_location is None:
        new_branch_pos_id, new_branch_pos = pick_new_branch_position(capped_surface)
    else:
        new_branch_pos_id, new_branch_pos = get_new_branch_position(branch_location, capped_surface)

    branch_to_manipulate_end_point = get_end_point(branch_to_manipulate)
    centerlines = update_centerlines(centerlines_complete, branch_to_manipulate_end_point)

    # Select diverging branch from Brancher and compare with diverging branch found manually
    centerlines_branched = vmtk_branch_extractor(centerlines_complete)
    diverging_centerline_branch = get_diverging_centerline_branch(centerlines_branched, branch_to_manipulate_end_point)

    # Clip Voronoi diagram into bend and remaining part of geometry
    print("-- Clipping Voronoi diagrams")
    voronoi_branch, voronoi_remaining = split_voronoi_with_centerlines(voronoi, [diverging_centerline_branch,
                                                                                 centerlines])

    # Get surface normals
    old_normal = get_estimated_surface_normal(diverging_centerline_branch)
    new_normal = get_exact_surface_normal(capped_surface, new_branch_pos_id)

    # Define rotation between surface normal vectors
    u, surface_normals_angle = get_rotation_axis_and_angle(new_normal, old_normal)
    R_u = get_rotation_matrix(u, surface_normals_angle)

    # Define rotation around new surface normal
    R_z = get_rotation_matrix(-new_normal, angle)

    # Define translation parameters
    dx, origo = get_translation_parameters(centerlines, diverging_centerline_branch, new_branch_pos)

    # Move branch centerline and voronoi diagram
    moved_and_rotated_voronoi_branch = move_and_rotate_voronoi_branch(voronoi_branch, dx, R_u, R_z, origo)
    moved_and_rotated_centerline_branch = move_and_rotate_centerline_branch(diverging_centerline_branch, origo, R_u,
                                                                            R_z, dx)
    # Create new voronoi diagram and new centerlines
    new_voronoi = merge_data([voronoi_remaining, moved_and_rotated_voronoi_branch])
    write_polydata(new_voronoi, new_voronoi_path)

    new_centerlines = merge_data([centerlines, moved_and_rotated_centerline_branch])
    write_polydata(new_centerlines, new_centerlines_path)

    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)
    write_polydata(new_surface, unprepared_output_filepath)

    print("-- Smoothing, clean, and check surface")
    new_surface = prepare_surface_output(new_surface, surface,
                                         new_centerlines, output_filepath,
                                         test_merge=False, changed=True,
                                         old_centerline=centerlines_complete)

    write_polydata(new_surface, output_filepath)


def check_branch_number(branch_to_manipulate_number, centerlines_complete):
    num_lines = centerlines_complete.GetNumberOfLines() - 1
    if branch_to_manipulate_number > num_lines:
        raise RuntimeError("\nERROR: Branch number cannot exceed number of centerlines." +
                           " Number of selectable centerlines for this model is {}.".format(num_lines))


def get_new_branch_position(branch_location, capped_surface):
    surface_locator = get_locator(capped_surface)
    new_branch_pos_id = surface_locator.FindClosestPoint(branch_location)
    new_branch_pos = capped_surface.GetPoint(new_branch_pos_id)

    return new_branch_pos_id, new_branch_pos


def pick_new_branch_position(capped_surface):
    print("\nPlease select new point to move branch in the render window.")
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(capped_surface)
    seed_selector.text = "Press space to select the point where you want to move the branch to, 'u' to undo.\n"
    seed_selector.Execute()
    new_branch_pos_id = seed_selector.GetTargetSeedIds().GetId(0)
    new_branch_pos = capped_surface.GetPoint(new_branch_pos_id)

    return new_branch_pos_id, new_branch_pos


def get_diverging_centerline_end(diverging_centerlines, diverging_id, region_points):
    print("-- Clipping a centerline divering in the region of interest.")
    patch_diverging_line = clip_diverging_line(extract_single_line(diverging_centerlines, 0),
                                               region_points[0], diverging_id)
    diverging_centerline = extract_single_line(patch_diverging_line, 1)
    diverging_centerline_end = get_end_point(diverging_centerline)

    return diverging_centerline_end


def get_end_point(centerline):
    centerline_end_point = centerline.GetPoint(centerline.GetNumberOfPoints() - 1)

    return centerline_end_point


def get_translation_parameters(centerlines, diverging_centerline_branch, new_branch_pos):
    locator = get_locator(centerlines)
    new_branch_pos_on_cl = centerlines.GetPoint(locator.FindClosestPoint(new_branch_pos))
    adjusted_branch_pos = (np.asarray(new_branch_pos) - np.asarray(new_branch_pos_on_cl)) * 0.8 + np.asarray(
        new_branch_pos_on_cl)
    dx = np.asarray(adjusted_branch_pos) - np.asarray(diverging_centerline_branch.GetPoint(0))
    origo = np.asarray(adjusted_branch_pos)

    return dx, origo


def get_exact_surface_normal(capped_surface, new_branch_pos_id):
    capped_surface_with_normals = vmtk_surface_normals(capped_surface)

    new_normal = capped_surface_with_normals.GetPointData().GetNormals().GetTuple(new_branch_pos_id)
    new_normal /= la.norm(new_normal)

    return new_normal


def get_diverging_centerline_branch(centerlines_branched, diverging_centerline_end):
    for i in range(centerlines_branched.GetNumberOfLines()):
        centerline_branch = extract_single_line(centerlines_branched, i)
        if get_end_point(centerline_branch) == diverging_centerline_end:
            return centerline_branch


def update_centerlines(centerlines, diverging_centerline_end):
    remaining_centerlines = []
    for i in range(centerlines.GetNumberOfLines()):
        diverging_centerline = extract_single_line(centerlines, i)
        if get_end_point(diverging_centerline) != diverging_centerline_end:
            remaining_centerlines.append(diverging_centerline)

    centerlines = merge_data(remaining_centerlines)

    return centerlines


def get_estimated_surface_normal(diverging_centerline_branch):
    line = vmtk_centerline_geometry(diverging_centerline_branch, True);
    curvature = get_array("Curvature", line)
    first_local_maxima_id = argrelextrema(curvature, np.greater)[0][0]

    # TODO: Generalize choice of end point factor
    factor = 0.6
    start_point = 0
    end_point = int(first_local_maxima_id * factor)
    normal_vector = np.asarray(diverging_centerline_branch.GetPoint(end_point)) - np.asarray(
        diverging_centerline_branch.GetPoint(start_point))

    return normal_vector / la.norm(normal_vector)


def move_and_rotate_centerline_branch(centerline_branch, origo, R_u, R_z, dx):
    # Locator to find closest point on centerline
    number_of_points = centerline_branch.GetNumberOfPoints()

    centerline = vtk.vtkPolyData()
    centerline_points = vtk.vtkPoints()
    centerline_cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, number_of_points)

    centerline_cell_array.InsertNextCell(number_of_points)
    radius_array_data = centerline_branch.GetPointData().GetArray(radiusArrayName).GetTuple1

    for p in range(centerline_branch.GetNumberOfPoints()):
        point = centerline_branch.GetPoint(p)

        # Translate
        point = np.asarray(point) + dx

        # Rotate
        point = np.dot(R_u, point - origo) + origo
        point = np.dot(R_z, point - origo) + origo

        centerline_points.InsertNextPoint(point)
        radius_array.SetTuple1(p, radius_array_data(p))
        centerline_cell_array.InsertCellPoint(p)

    centerline.SetPoints(centerline_points)
    centerline.SetLines(centerline_cell_array)
    centerline.GetPointData().AddArray(radius_array)

    return centerline


def get_rotation_axis_and_angle(n_new, n_old):
    z = np.asarray(n_old)
    z_prime = np.asarray(n_new)

    u = np.cross(z, z_prime)
    u /= la.norm(u)

    angle = np.arccos(z.dot(z_prime) / (la.norm(z) * la.norm(z_prime)))

    return u, angle


def get_rotation_matrix(u, angle):
    # Set up rotation matrix using Euler–Rodrigues formula
    u_cross_matrix = np.asarray([[0, -u[2], u[1]],
                                 [u[2], 0, -u[0]],
                                 [-u[1], u[0], 0]])
    u_outer = np.outer(u, u)
    R = np.cos(angle) * np.eye(3) + np.sin(angle) * u_cross_matrix + (1 - np.cos(angle)) * u_outer

    return R


def move_and_rotate_voronoi_branch(voronoi, dx, R_u, R_z, origo):
    # Voronoi diagram
    n = voronoi.GetNumberOfPoints()
    new_voronoi = vtk.vtkPolyData()
    voronoi_points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, n)

    # Iterate through Voronoi diagram and manipulate
    for i in range(n):
        point = voronoi.GetPoint(i)
        point_radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

        # Translate and rotate
        point = np.asarray(point) + dx
        point = np.dot(R_u, point - origo) + origo
        point = np.dot(R_z, point - origo) + origo

        # Change radius
        radius_array.SetTuple1(i, point_radius)
        voronoi_points.InsertNextPoint(point)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i)

    new_voronoi.SetPoints(voronoi_points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)

    return new_voronoi


def get_sorted_lines(centerlines_complete):
    lines = [extract_single_line(centerlines_complete, i) for i in range(centerlines_complete.GetNumberOfLines())]
    sorted_lines = sorted(lines, key=functools.cmp_to_key(compare_to))

    return sorted_lines


def get_branch(branch_to_manipulate_number, centerlines_complete, longest_centerline):
    break_points = dict()
    for i in range(centerlines_complete.GetNumberOfLines()):
        current_line = extract_single_line(centerlines_complete, i)
        tolerance = get_tolerance(current_line)
        for j in range(current_line.GetNumberOfPoints()):
            p1 = np.asarray(current_line.GetPoint(j))
            p2 = np.asarray(longest_centerline.GetPoint(j))
            if distance(p1, p2) > tolerance:
                break_points[j] = i
                break

    break_points_keys = sorted(break_points.keys())
    selected_branch_id = break_points[break_points_keys[branch_to_manipulate_number - 1]]
    selected_branch = extract_single_line(centerlines_complete, selected_branch_id)

    return selected_branch


def compare_to(line0, line1):
    len0 = len(get_curvilinear_coordinate(line0))
    len1 = len(get_curvilinear_coordinate(line1))
    if len0 > len1:
        return 1
    return -1


def pick_branch(capped_surface, centerlines):
    print("\nPlease select the branch you with to move in the render window.")
    # Select point on surface
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(capped_surface)
    seed_selector.text = "Press space to select the branch you want to move, 'u' to undo.\n"
    seed_selector.Execute()
    branch_surface_point_id = seed_selector.GetTargetSeedIds().GetId(0)
    surface_point = capped_surface.GetPoint(branch_surface_point_id)

    # Find closest centerline
    closest_dist = 1e9
    branch_to_manipulate = None
    for line_to_compare in centerlines:
        locator = get_locator(line_to_compare)
        closest_point_id = locator.FindClosestPoint(surface_point)
        closest_point = line_to_compare.GetPoint(closest_point_id)
        dist = distance(closest_point, surface_point)
        if dist < closest_dist:
            closest_dist = dist
            branch_to_manipulate = line_to_compare

    return branch_to_manipulate


def read_command_line():
    """
    Read arguments from commandline
    """
    # Description of the script
    description = "Moves a selected part of a tubular geometry, " + \
                  "in two (horizontal, vertical) geometry-specific directions. " + \
                  "Magnitude of movement is defined by the parameters alpha and beta" + \
                  "Primary script used for application in blood vessels."

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    # Add common arguments
    add_common_arguments(parser)

    parser.add_argument('-bl', "--branch-loc", nargs="+", type=float, default=None, metavar="branch_location",
                        help="If this parameter is provided, the branch to be manipulated will be moved to the point "
                             "on the surface closest to this point. Example providing the point (1, 5, -1):" +
                             " --region-points 1 5 -1")

    # Arguments for rotation
    parser.add_argument('-a', '--angle', type=float, default=0,
                        help="The manipulated branch is rotated an angle 'a' around the new" +
                             " surface normal vector. 'a' is assumed to be in degrees," +
                             " and not radians. Default is no rotation.", metavar="surface_normal_axis_angle")
    # Argument for selecting branch
    parser.add_argument('-bn', '--branch-number', type=int, default=0,
                        help="The number corresponding the branch to manipulate. " +
                             "The branches are ordered from 1 to N, " +
                             "from upstream to downstream, relative to the inlet. " +
                             "If not selected, the user must manually select the branch "
                             "to manipulate. ", metavar="branch_number")

    # Parse paths to get default values
    args = parser.parse_args()
    angle_to_radians = args.angle * math.pi / 180  # Convert from deg to rad

    # Output file argument
    args = parser.parse_args()

    if args.no_smooth_point is not None and len(args.no_smooth_point):
        if len(args.no_smooth_point) % 3 != 0:
            raise ValueError("ERROR: Please provide the no smooth point(s) as a multiple of 3")

    return dict(input_filepath=args.ifile, smooth=args.smooth, smooth_factor=args.smooth_factor,
                output_filepath=args.ofile, poly_ball_size=args.poly_ball_size,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point,
                resampling_step=args.resampling_step, angle=angle_to_radians,
                branch_to_manipulate_number=args.branch_number, branch_location=args.branch_loc)


def main_branch():
    manipulate_branch(**read_command_line())


if __name__ == "__main__":
    manipulate_branch(**read_command_line())

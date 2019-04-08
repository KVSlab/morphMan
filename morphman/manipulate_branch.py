##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from common import *


def test_data(capped_surface):
    branch_point_surface_id = 12728
    region_points = [[37.119407653808594, 27.06595230102539, 33.817745208740234],
                     [37.82146072387695, 22.682409286499023, 28.8624267578125]]
    branch_point = capped_surface.GetPoint(branch_point_surface_id)
    return branch_point_surface_id, region_points, branch_point


def get_new_branch_position(capped_surface):
    print("\nPlease select point to move branch in the render window.")
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(capped_surface)
    seed_selector.text = "Press space to select the point where you want to move the branch to, 'u' to undo.\n"
    seed_selector.Execute()
    new_branch_pos_id = seed_selector.GetTargetSeedIds().GetId(0)
    new_branch_pos = capped_surface.GetPoint(new_branch_pos_id)

    return new_branch_pos_id, new_branch_pos


def manipulate_branch(input_filepath, output_filepath, smooth, smooth_factor, region_of_interest, region_points,
                      poly_ball_size, no_smooth, no_smooth_point, resampling_step, angle):
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
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'landmarking']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        poly_ball_size (list): Resolution of polyballs used to create surface.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
        resampling_step (float): Resampling length when resampling centerline.
        angle (float): Angle to rotate around new surface normal vector. Default i no rotation.
    """

    # Input filenames
    base_path = get_path_names(input_filepath)

    # Centerlines filenames
    centerlines_path = base_path + "_centerline.vtp"
    diverging_centerline_branch_path = base_path + "_centerlines_branch.vtp"
    new_centerlines_path = base_path + "_centerline_moved_and_rotated.vtp"

    # Voronoi diagrams filenames
    new_voronoi_path = base_path + "_voronoi_moved_and_rotated.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Compute centerlines
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

    # TODO: For test only
    new_branch_pos_id, region_points, new_branch_pos = test_data(capped_surface)

    if region_points is None:
        # Get region of interest
        _, _, _, region_points, _ = get_line_to_change(capped_surface, centerlines_complete,
                                                       region_of_interest, "bend", region_points, 0.0)
        region_points = [[region_points[3 * i], region_points[3 * i + 1], region_points[3 * i + 2]]
                         for i in range(len(region_points) // 3)]

        # Get new position of branch on model surface
        new_branch_pos_id, new_branch_pos = get_new_branch_position(capped_surface)

    # Branch-extractor
    Brancher = vmtkscripts.vmtkBranchExtractor()
    Brancher.Centerlines = centerlines_complete
    Brancher.RadiusArrayName = radiusArrayName
    Brancher.Execute()
    centerlines_branched = Brancher.Centerlines

    # Set and get clipping points, centerlines and diverging centerlines
    centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        find_region_of_interest_and_diverging_centerlines(centerlines_complete, region_points)
    diverging_centerline_ispresent = diverging_centerlines is not None
    diverging_id = None if len(diverging_ids) == 0 else diverging_ids[0]

    # Handle diverging centerlines within region of interest
    if diverging_centerline_ispresent:
        print("-- Clipping a centerline divering in the region of interest.")
        patch_diverging_line = clip_diverging_line(extract_single_line(diverging_centerlines, 0),
                                                   region_points[0], diverging_id)
        diverging_centerline = extract_single_line(patch_diverging_line, 1)
        diverging_centerline_end = diverging_centerline.GetPoint(diverging_centerline.GetNumberOfPoints() - 1)
    else:
        raise RuntimeError("No diverging branch detected! Cannot translate nothing.")

    # Select diverging branch from Brancher and compare with diverging branch found manually
    for i in range(centerlines_branched.GetNumberOfLines()):
        centerline_branch = extract_single_line(centerlines_branched, i)
        if centerline_branch.GetPoint(centerline_branch.GetNumberOfPoints() - 1) == diverging_centerline_end:
            diverging_centerline_branch = centerline_branch
            write_polydata(diverging_centerline_branch, diverging_centerline_branch_path)
            break

    diverging_centerline_branch_start = diverging_centerline_branch.GetPoint(0)
    radiusArrayDiv = diverging_centerline_branch.GetPointData().GetArray(radiusArrayName)
    on_div_cl_misr = radiusArrayDiv.GetTuple1(0)
    end_point, rad, _ = move_past_sphere(diverging_centerline_branch, diverging_centerline_branch_start, on_div_cl_misr,
                                         start=0, stop=diverging_centerline_branch.GetNumberOfPoints(), step=1, X=0.5)

    # Clip Voronoi diagram into bend and remaining part of geometry
    print("-- Clipping Voronoi diagrams")
    voronoi_branch, voronoi_remaining = split_voronoi_with_centerlines(voronoi,
                                                                       [diverging_centerline_branch,
                                                                        centerlines])
    # Get surface normals
    SurfaceNormals = vmtkscripts.vmtkSurfaceNormals()
    SurfaceNormals.Surface = capped_surface
    SurfaceNormals.NormalsArrayName = surfaceNormalsArrayName

    SurfaceNormals.Execute()
    capped_surface_with_normals = SurfaceNormals.Surface
    new_normal = capped_surface_with_normals.GetPointData().GetNormals().GetTuple(new_branch_pos_id)
    # new_normal = capped_surface_with_normals.GetPointData().GetNormals().GetTuple(branch_point_id.GetId(0))
    new_normal /= la.norm(new_normal)

    old_normal = get_estimated_normal_vector(diverging_centerline_branch)
    old_normal /= la.norm(old_normal)

    # Define rotation between surface normal vectors
    u, surface_normals_angle = get_rotation_axis_and_angle(new_normal, old_normal)
    R_u = get_rotation_matrix(u, surface_normals_angle)

    # Define rotation around new surface normal
    R_z = get_rotation_matrix(-new_normal, angle)

    # Define translation parameters
    dx = np.asarray(new_branch_pos) - np.asarray(diverging_centerline_branch.GetPoint(0))
    origo = np.asarray(new_branch_pos)

    # Move branch centerline and voronoi diagram
    moved_voronoi_branch = move_voronoi_branch(voronoi_branch, dx)
    rotated_voronoi_branch_0 = rotate_voronoi_branch(moved_voronoi_branch, origo, R_u)
    rotated_voronoi_branch = rotate_voronoi_branch(rotated_voronoi_branch_0, origo, R_z)

    moved_and_rotated_centerline_branch = move_and_rotate_centerline_branch(diverging_centerline_branch, origo, R_u,
                                                                            R_z, dx)

    # Create new voronoi diagram and new centerlines
    new_voronoi = merge_data([voronoi_remaining, rotated_voronoi_branch])
    write_polydata(new_voronoi, new_voronoi_path)

    new_centerlines = merge_data([centerlines, moved_and_rotated_centerline_branch])
    write_polydata(new_centerlines, new_centerlines_path)

    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    print("-- Smoothing, clean, and check surface")
    new_surface = prepare_surface_output(new_surface, surface,
                                         new_centerlines, output_filepath,
                                         test_merge=True, changed=True,
                                         old_centerline=centerlines_complete)

    write_polydata(new_surface, output_filepath)


def get_estimated_normal_vector(diverging_centerline_branch):
    line = vmtk_centerline_geometry(diverging_centerline_branch, True);
    curvature = get_array("Curvature", line)
    first_local_maxima_id = argrelextrema(curvature, np.greater)[0][0]

    # TODO: Generalize choice of end point factor
    factor = 0.6
    start_point = 0
    end_point = int(first_local_maxima_id * factor)
    normal_vector = np.asarray(diverging_centerline_branch.GetPoint(end_point)) - np.asarray(
        diverging_centerline_branch.GetPoint(start_point))

    return normal_vector


def move_and_rotate_centerline_branch(centerline_branch, origo, R_u, R_z, dx):
    origo = np.asarray(origo)

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


def rotate_voronoi_branch(voronoi, origo, R):
    origo = np.asarray(origo)

    # Voronoi diagram
    n = voronoi.GetNumberOfPoints()
    new_voronoi = vtk.vtkPolyData()
    voronoi_points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, n)

    # Iterate through Voronoi diagram and manipulate
    k = 0
    for i in range(n):
        point = voronoi.GetPoint(i)
        point_radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

        # Rotate
        point = np.dot(R, point - origo) + origo

        # Change radius
        radius_array.SetTuple1(k, point_radius)
        voronoi_points.InsertNextPoint(point)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(k)
        k += 1

    new_voronoi.SetPoints(voronoi_points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)
    return new_voronoi


def move_voronoi_branch(voronoi, dx):
    # Voronoi diagram
    n = voronoi.GetNumberOfPoints()
    new_voronoi = vtk.vtkPolyData()
    voronoi_points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, n)

    # Iterate through Voronoi diagram and manipulate
    k = 0
    for i in range(n):
        point = voronoi.GetPoint(i)
        point_radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

        # Translate
        point = np.asarray(point) + dx

        # Change radius
        radius_array.SetTuple1(k, point_radius)
        voronoi_points.InsertNextPoint(point)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(k)
        k += 1

    new_voronoi.SetPoints(voronoi_points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)
    return new_voronoi


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

    # Set region of interest:
    parser.add_argument("-r", "--region-of-interest", type=str, default="manual",
                        choices=["manual", "commandline", "landmarking"],
                        help="The method for defining the region to be changed. There are" +
                             " three options: 'manual', 'commandline', 'landmarking'. In" +
                             " 'manual' the user will be provided with a visualization of the" +
                             " input surface, and asked to provide an end and start point of the" +
                             " region of interest. Note that not all algorithms are robust over" +
                             " bifurcations. If 'commandline' is provided, then '--region-points'" +
                             " is expected to be provided. Finally, if 'landmarking' is" +
                             " given, it will look for the output from running automated_landmarking.py.")

    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" +
                             " argument have to be given. The method expects two points" +
                             " which defines the start and end of the region of interest. " +
                             " Example providing the points (1, 5, -1) and (2, -4, 3):" +
                             " --region-points 1 5 -1 2 -4 3")

    # Arguments for rotation
    parser.add_argument('-a', '--angle', type=float, default=0,
                        help="The manipulated branch is rotated an angle 'a' around the new" +
                             " surface normal vector. 'a' is assumed to be in degrees," +
                             " and not radians. Default is no rotation.", metavar="surface_normal_axis_angle")

    # Parse paths to get default values
    args = parser.parse_args()
    angle_to_radians = args.angle * math.pi / 180  # Convert from deg to rad

    # Output file argument
    args = parser.parse_args()

    if args.no_smooth_point is not None and len(args.no_smooth_point):
        if len(args.no_smooth_point) % 3 != 0:
            raise ValueError("ERROR: Please provide the no smooth point(s) as a multiple of 3")

    return dict(input_filepath=args.ifile, smooth=args.smooth,
                smooth_factor=args.smooth_factor,
                output_filepath=args.ofile, poly_ball_size=args.poly_ball_size,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point,
                region_of_interest=args.region_of_interest, region_points=args.region_points,
                resampling_step=args.resampling_step, angle=angle_to_radians)


def main_branch():
    manipulate_branch(**read_command_line())


if __name__ == "__main__":
    manipulate_branch(**read_command_line())

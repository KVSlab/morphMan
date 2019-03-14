##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from common import *

# Local import

surfaceNormalsArrayName = 'SurfaceNormalArray'
frenetTangentArrayName = 'FrenetTangent'


def manipulate_branch(input_filepath, output_filepath, smooth, smooth_factor, region_of_interest, region_points,
                      poly_ball_size, no_smooth, no_smooth_point, resampling_step):
    """
    Primary script for moving a selected part of any blood vessel.
    Relies on an input centerline, a surface geometry of a 3D blood vessel network,
    and the magnitude of movement (alpha, beta) in horizontal and vertical direction respectively.

    Defines in- and out-put files, computes voronoi diagram, and
    moves voronoi diagram in the horizontal direction based on the
    definitions in "Investigating the Interaction Between Morphology
    of the Anterior Bend and Aneurysm Initiation" (2018).
    Proceeds by computing the corresponding centerline in the new
    geometries, used for postprocessing, geometric analysis and
    meshing.

    Continuation in the move_vessel_vertically-method for vertical movement.

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
    """

    # Input filenames
    base_path = get_path_names(input_filepath)

    # Centerlines filenames
    centerlines_path = base_path + "_centerline.vtp"
    diverging_centerline_branch_path = base_path + "_centerlines_branch.vtp"

    # Voronoi diagrams filenames
    new_voronoi_path = base_path + "_voronoi_interpolated.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Compute centerlines
    inlet, outlets = get_centers(surface, base_path)

    print("-- Compute centerlines and Voronoi diagram")
    centerlines_complete, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                                  capped_surface, resampling=resampling_step,
                                                                  smooth=False, base_path=base_path)
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines_complete, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)

        # FIXME: For testing only:
        # write_polydata(centerlines_complete, "COMPLETE.vtp")
        # branch_point = (52.46571350097656, 28.395702362060547, 17.509746551513672)

    # region_points = np.asarray([[34.54340362548828, 7.8101677894592285, 43.735019683837899],
    #                            [36.73946762084961, 28.79076385498047, 36.41799545288086]])
    region_points = np.asarray([[38.820396423339844, 20.89110565185547, 28.940338134765625],
                                [38.287437438964844, 27.811443328857422, 34.812583923339844]])
    branch_point_surface_id = 6309

    if region_points is None:
        # Get region of interest
        _, _, _, region_points, _ = get_line_to_change(capped_surface, centerlines_complete,
                                                       region_of_interest, "bend", region_points, 0.0)
        region_points = [[region_points[3 * i], region_points[3 * i + 1], region_points[3 * i + 2]]
                         for i in range(len(region_points) // 3)]
        # Get place to put branch
        print("\nPlease select point to move branch in the render window.")
        seed_selector = vmtkPickPointSeedSelector()
        seed_selector.SetSurface(capped_surface)
        seed_selector.text = "Press space to select the point where you want to move the branch to, 'u' to undo.\n"
        seed_selector.Execute()
        branch_point_id = seed_selector.GetTargetSeedIds()
        branch_point = capped_surface.GetPoint(branch_point_id.GetId(0))
        branch_point_surface_id = branch_point_id.GetId(0)
        print(branch_point_surface_id)

    # Test branch-extractor
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
    diverging_point = centerlines.GetPoint(diverging_id)

    # Map branch surface point to closest centerline point
    locator = get_locator(centerlines)
    branch_point_on_cl_id = locator.FindClosestPoint(branch_point)
    branch_point_on_cl = centerlines.GetPoint(branch_point_on_cl_id)

    # Handle diverging centerlines within region of interest
    if diverging_centerline_ispresent:
        print("-- Clipping a centerline divering in the region of interest.")
        patch_diverging_line = clip_diverging_line(extract_single_line(diverging_centerlines, 0),
                                                   region_points[0], diverging_id)
        diverging_centerline = extract_single_line(patch_diverging_line, 1)
        diverging_centerline_end = diverging_centerline.GetPoint(diverging_centerline.GetNumberOfPoints() - 1)
    else:
        raise RuntimeError("No diverging branch detected! Cannot translate nothing.")

    # Select diverging branch
    for i in range(centerlines_branched.GetNumberOfLines()):
        centerline_branch = extract_single_line(centerlines_branched, i)
        if centerline_branch.GetPoint(centerline_branch.GetNumberOfPoints() - 1) == diverging_centerline_end:
            diverging_centerline_branch = centerline_branch
            write_polydata(diverging_centerline_branch, diverging_centerline_branch_path)
            break

    locator = get_locator(centerlines)
    diverging_centerline_branch_start = diverging_centerline_branch.GetPoint(0)
    on_cl_id = locator.FindClosestPoint(diverging_centerline_branch_start)
    on_cl_point = centerlines.GetPoint(on_cl_id)
    radiusArray = centerlines.GetPointData().GetArray(radiusArrayName)

    locator = get_locator(diverging_centerline_branch)
    radiusArrayDiv = diverging_centerline_branch.GetPointData().GetArray(radiusArrayName)
    on_div_cl_misr = radiusArrayDiv.GetTuple1(0)
    end_point, rad, _ = move_past_sphere(diverging_centerline_branch, diverging_centerline_branch_start, on_div_cl_misr,
                                         start=0, stop=diverging_centerline_branch.GetNumberOfPoints(), step=1, X=0.5)
    threshold_misr = distance(end_point, on_cl_point)
    on_div_cl_id = locator.FindClosestPoint(end_point)

    # Clip Voronoi diagram into bend and remaining part of geometry
    print("-- Clipping Voronoi diagrams")
    voronoi_branch, voronoi_remaining = split_voronoi_with_centerlines(voronoi,
                                                                       [diverging_centerline_branch,
                                                                        centerlines])

    # Get surface normals
    SurfaceNormals = vmtkscripts.vmtkSurfaceNormals()
    SurfaceNormals.Surface = capped_surface
    # SurfaceNormals.NormalsArrayName = surfaceNormalsArrayName
    SurfaceNormals.NormalsArrayName = "SurfaceNormalArray"
    SurfaceNormals.Execute()
    capped_surface_with_normals = SurfaceNormals.Surface
    normal = capped_surface_with_normals.GetPointData().GetNormals().GetTuple(branch_point_surface_id)
    # normal = capped_surface_with_normals.GetPointData().GetNormals().GetTuple(branch_point_id.GetId(0))
    normal = normal / la.norm(normal)

    # TODO: TESTING
    testing = False
    if testing:
        surface_no_branch = create_new_surface(voronoi_remaining, poly_ball_size=[250, 250, 250])
        SurfaceNormals = vmtkscripts.vmtkSurfaceNormals()
        SurfaceNormals.Surface = surface_no_branch
        SurfaceNormals.NormalsArrayName = surfaceNormalsArrayName
        SurfaceNormals.Execute()
        capped_surface_with_normals = SurfaceNormals.Surface

        surface_locator = get_locator(surface_no_branch)
        surface_id = surface_locator.FindClosestPoint(diverging_centerline_branch.GetPoint(0))
        normal_old = capped_surface_with_normals.GetPointData().GetNormals().GetTuple(surface_id)
        normal_old = normal_old / la.norm(normal_old)
    else:
        normal_old = np.asarray(diverging_centerline_branch.GetPoint(5)) - np.asarray(
            diverging_centerline_branch.GetPoint(0))
        normal_old = normal_old / la.norm(normal)

    # Define translation parameters
    locator = get_locator(centerlines)
    branch_point_on_cl_id_old = locator.FindClosestPoint(diverging_centerline_branch.GetPoint(0))
    R_u, R_z = get_rotation_matrices(normal, normal_old, branch_point_on_cl_id_old, branch_point_on_cl_id,
                                     centerlines_complete)
    diverging_centerline_branch_clip = extract_single_line(diverging_centerline_branch, 0, startID=on_div_cl_id)
    dx = np.asarray(branch_point) - np.asarray(diverging_centerline_branch_clip.GetPoint(0))
    origo = np.asarray(branch_point)

    # Move branch centerline and voronoi diagram
    moved_voronoi_branch = move_voronoi_branch(voronoi_branch, centerlines, threshold_misr, origo, R_u, R_z, dx)
    centerline_branch_moved = move_centerline_branch(diverging_centerline_branch_clip, origo, R_u, R_z, dx)
    write_polydata(moved_voronoi_branch, "ABRANCH.vtp")
    write_polydata(centerline_branch_moved, "ALINE.vtp")

    file = open("normal.particles", "w")

    for i in range(100):
        p = np.asarray(branch_point) + i * normal
        file.write("%s %s %s\n" % (p[0], p[1], p[2]))
    file.close()

    file = open("old_normal.particles", "w")

    for i in range(100):
        p = diverging_centerline_branch.GetPoint(0) + i * normal_old
        file.write("%s %s %s\n" % (p[0], p[1], p[2]))
    file.close()

    # Set interpolation starting point
    locator = get_locator(centerlines)
    moved_centerline_branch_start = centerline_branch_moved.GetPoint(0)
    moved_branch_point_cl_id = locator.FindClosestPoint(moved_centerline_branch_start)
    moved_branch_point_cl_misr = radiusArray.GetTuple1(moved_branch_point_cl_id)
    moved_branch_point_cl_point = centerlines.GetPoint(moved_branch_point_cl_id)
    cl_start_point, start_rad, _ = move_past_sphere(centerlines, moved_branch_point_cl_point,
                                                    moved_branch_point_cl_misr,
                                                    start=moved_branch_point_cl_id, stop=0, step=-1, X=3.0)

    # Set additional point for interpolation
    p_div = np.asarray(moved_centerline_branch_start)
    p_start = np.asarray(cl_start_point)
    p_mid = np.asarray(moved_branch_point_cl_point)
    interp_point_start = 1 / 3 * p_start + 1. / 3 * p_div + 1 / 3 * p_mid

    # Clip old centerline to define up- and down-stream centerlines relative to new branch
    clipping_points = vtk.vtkPoints()
    for p in [p_start, p_start]:
        clipping_points.InsertNextPoint(p)
    patch_cl_complete = create_parent_artery_patches(centerlines, clipping_points, siphon=True)

    patch_cl_upstream = extract_single_line(patch_cl_complete, 0)
    patch_cl_downstream = merge_data(
        [extract_single_line(patch_cl_complete, i) for i in range(1, patch_cl_complete.GetNumberOfLines())])

    patch_upstream = merge_data([patch_cl_upstream, centerline_branch_moved])

    # Interpolate centerline between upstream centerline and moved branch
    print("-- Interpolate centerlines.")
    interpolated_cl_upstream = interpolate_patch_centerlines(patch_upstream, diverging_centerlines,
                                                             None,
                                                             None, False)
    # Split remaining Voronoi diagram relative to branch
    voronoi_start, voronoi_end = split_voronoi_with_centerlines(voronoi_remaining,
                                                                [patch_cl_upstream,
                                                                 patch_cl_downstream])

    write_polydata(interpolated_cl_upstream, "interpolated_centerline_to_movespot.vtp")
    # Interpolate voronoi diagram
    print("-- Interpolate voronoi diagram.")
    p1 = np.asarray(clipping_points.GetPoint(0))
    p2 = np.asarray(moved_centerline_branch_start)
    vtk_points = vtk.vtkPoints()
    for p in [p1, p2]:
        vtk_points.InsertNextPoint(p)
    clipping_points = (vtk_points, np.asarray([p1, p2]))
    patch_cl_upstream_and_branch = merge_data([patch_cl_upstream, centerline_branch_moved])
    interpolated_voronoi = interpolate_voronoi_diagram(interpolated_cl_upstream, patch_cl_upstream_and_branch,
                                                       merge_data([voronoi_start, moved_voronoi_branch]),
                                                       clipping_points,
                                                       [], cylinder_factor=7.0)
    write_polydata(interpolated_voronoi, "voroINTER.vtp")
    # Write a new surface from the new voronoi diagram
    print("-- Create new surface.")
    # new_voronoi = merge_data([voronoi_end, merge_data([voronoi_start, moved_voronoi_branch])])
    new_voronoi = merge_data([voronoi_end, interpolated_voronoi])
    write_polydata(new_voronoi, new_voronoi_path)

    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)
    write_polydata(new_surface, output_filepath)

    # print("-- Smoothing, clean, and check surface")
    """
    new_centerlines = merge_data([interpolated_cl_upstream, centerlines])
    new_surface = prepare_surface_output(new_surface, surface,
                                         new_centerlines, output_filepath,
                                         test_merge=True, changed=True,
                                         old_centerline=merge_data([centerlines,
                                                                    diverging_centerlines]))
    """


def reverse_centerline(centerline):
    number_of_points = centerline.GetNumberOfPoints()

    centerline_reverse = vtk.vtkPolyData()
    centerline_points = vtk.vtkPoints()
    centerline_cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, number_of_points)

    centerline_cell_array.InsertNextCell(number_of_points)
    radius_array_data = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1
    i = 0
    for p in range(number_of_points - 1, -1, -1):
        point = centerline.GetPoint(p)
        centerline_points.InsertNextPoint(np.asarray(point))
        radius_array.SetTuple1(i, radius_array_data(p))
        centerline_cell_array.InsertCellPoint(i)
        i += 1

    centerline_reverse.SetPoints(centerline_points)
    centerline_reverse.SetLines(centerline_cell_array)
    centerline_reverse.GetPointData().AddArray(radius_array)

    return centerline_reverse


def move_centerline_branch(centerline_branch, origo, R_u, R_z, dx):
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


def get_rotation_matrices(n, n_old, id, id_old, centerlines):
    z_prime = np.asarray(n)
    z = np.asarray(n_old)

    # Find tangents along centerline
    centerlines = vmtk_centerline_geometry(centerlines, False)
    frenet_tangent_array = get_array(frenetTangentArrayName, centerlines, k=3)
    t_0 = frenet_tangent_array[id_old]
    t_1 = frenet_tangent_array[id]

    # Define two coordinate systems
    y = np.cross(z, t_0)
    y_prime = np.cross(z_prime, t_1)
    x = np.cross(y, z)
    x_prime = np.cross(y_prime, z_prime)

    u = np.cross(z, z_prime)
    u = u / la.norm(u)

    u_cross_matrix = np.asarray([[0, -u[2], u[1]],
                                 [u[2], 0, -u[0]],
                                 [-u[1], u[0], 0]])

    angle = np.arccos(z.dot(z_prime) / (la.norm(z) * la.norm(z_prime)))

    u_outer = np.outer(u, u)

    # Set up rotation matrices
    R_u = np.cos(angle) * np.eye(3) + np.sin(angle) * u_cross_matrix + (1 - np.cos(angle)) * u_outer
    x_new = R_u.dot(x)

    theta = np.arccos(x_new.dot(x_prime) / (la.norm(x_new) * la.norm(x_prime)))
    cos_t = np.cos(theta + np.pi)
    sin_t = np.sin(theta + np.pi)

    R_z = np.asarray([[cos_t, -sin_t, 0],
                      [sin_t, cos_t, 0],
                      [0, 0, 1]])
    return R_u, R_z


def move_voronoi_branch(voronoi, centerline, threshold_misr, origo, R_u, R_z, dx):
    origo = np.asarray(origo)

    # Iterate through Voronoi diagram and manipulate

    # Voronoi diagram
    n = voronoi.GetNumberOfPoints()
    new_voronoi = vtk.vtkPolyData()
    voronoi_points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, n)

    # Iterate through Voronoi diagram and manipulate
    k = 0
    locator = get_locator(centerline)
    for i in range(n):
        point = voronoi.GetPoint(i)
        closest_cl_point = centerline.GetPoint(locator.FindClosestPoint(point))
        dist = distance(point, closest_cl_point)
        if dist < threshold_misr:
            continue

        point_radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

        # Translate
        point = np.asarray(point) + dx

        # Rotate
        point = np.dot(R_u, point - origo) + origo
        point = np.dot(R_z, point - origo) + origo

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
    # "Variation" arguments
    # parser.add_argument("--alpha", type=float, default=0.0,
    #                    help="Compression factor in vertical direction, " +
    #                         "ranging from -1.0 to 1.0, defining the magnitude " +
    #                         "of stretching or compression of the tubular structure.")
    # parser.add_argument("--beta", type=float, default=0.0,
    #                    help="Compression factor in vertical direction,  " +
    #                         "ranging from -1.0 to 1.0, defining the magnitude " +
    #                         "of stretching or compression of the tubular structure.")
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
                resampling_step=args.resampling_step)


def main_branch():
    manipulate_branch(**read_command_line())


if __name__ == "__main__":
    manipulate_branch(**read_command_line())

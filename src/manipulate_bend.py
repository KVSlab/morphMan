from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from common import *


def move_vessel(input_filepath, output_filepath, smooth, smooth_factor, region_of_interest, region_points,
                alpha, beta, poly_ball_size, no_smooth, no_smooth_point):
    """
    Primary script for moving a selected part of any blood vessel.
    Relies on an input centerline, a surface geometry of a 3D blood vessel network,
    and the magnitude of movement (alpha, beta) in horizontal and vertical direction respectively.

    Defines in- and out-put files, computes voronoi diagram, and
    moves voronoi diagram in the horizontal direction based on the
    definitions in "Investigating the Interaction Between Morphology
    of the Anterior Bend and Aneurysm Initiation" (2018).
    Procedes by computing the corresponding centerline in the new
    geometries, used for postprocessing, geometric analysis and
    meshing.

    Continuation in the move_vessel_vertically-method for vertical movement.

    Args:
        input_filepath (str): Path to input surface.
        output_filepath (str): Path to output the manipulated surface.
        smooth (bool): Smooth the Voronoi diagram.
        smooth_factor (float): Smoothing factor used for Voronoi diagram smoothing.
        region_of_interest (str): Method for setting the region of interest ['manuall' | 'commandline' | 'landmarking']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        poly_ball_size (list): Resolution of polyballs used to create surface.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
    """

    # Input filenames
    base_path = get_path_names(input_filepath)

    # Centerlines
    centerline_clipped_path = base_path + "_centerline_clipped.vtp"
    centerline_clipped_part_path = base_path + "_centerline_clipped_part.vtp"
    new_centerlines_path = base_path + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines_path_tmp = base_path + "new_centerlines_alpha_%s_beta_%s_tmp.vtp" % (alpha, beta)
    centerlines_path = base_path + "_centerline.vtp"

    # Voronoi diagrams
    voronoi_remaining_path = base_path + "_voronoi_remaining.vtp"
    voronoi_siphon_path = base_path + "_voronoi_siphon.vtp"

    # Surface information
    point_path = base_path + "_carotid_siphon_points.particles"  # Hard coded for consistency

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Compute centerlines
    inlet, outlets = get_centers(surface, base_path)

    print("Compute centerlines and Voronoi diagram")
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=0.1,
                                                         smooth=False, base_path=base_path)
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)

    # Get region of interest
    if region_of_interest == "landmarking":
        if not path.exists(point_path):
            RuntimeError(("The given .particles file: %s does not exist. Please run" + \
                          " landmarking with automated_geometric_quantities.py first.") % point_path)
        region_points = np.loadtxt(point_path)
    else:
        _, region_points = get_line_to_change(capped_surface, centerlines,
                                              region_of_interest, "bend", region_points, 0)
        region_points = [[region_points[3 * i], region_points[3 * i + 1], region_points[3 * i + 2]]
                         for i in range(len(region_points) // 3)]

    # Set and get clipping points, centerlines and diverging centerlines
    centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        find_region_of_interest_and_diverging_centerlines(centerlines, region_points)
    diverging_centerline_ispresent = diverging_centerlines is not None
    diverging_id = None if len(diverging_ids) == 0 else diverging_ids[0]

    # Handle diverging centerlines within region of interest
    if diverging_centerline_ispresent:
        print("Clipping opthamlic artery")
        patch_diverging_line = clip_diverging_line(extract_single_line(diverging_centerlines, 0), region_points[0],
                                                   diverging_id)

    # Clip centerline
    print("Clipping centerlines.")
    locator = get_locator(extract_single_line(centerlines, 0))
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    p1 = centerlines.GetPoint(id1)
    p2 = centerlines.GetPoint(id2)
    centerline_remaining = CreateParentArteryPatches(centerlines,
                                                     region_points_vtk, siphon=True)
    centerline_siphon = extract_single_line(centerlines, 0, startID=id1, endID=id2)

    if diverging_centerline_ispresent:
        eyeline_end = extract_single_line(patch_diverging_line, 1)
        centerline_siphon = merge_data([centerline_siphon, eyeline_end])

    write_polydata(centerline_remaining, centerline_clipped_path)
    write_polydata(centerline_siphon, centerline_clipped_part_path)

    # Clip Voronoi diagram into
    # siphon and remaining part of geometry
    print("Clipping Voronoi diagrams")
    if not path.exists(voronoi_siphon_path):
        masked_voronoi_clip = MaskVoronoiDiagram(voronoi, centerline_siphon)
        voronoi_siphon = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi_clip)
        write_polydata(voronoi_siphon, voronoi_siphon_path)
    else:
        voronoi_siphon = read_polydata(voronoi_siphon_path)

    if not path.exists(voronoi_remaining_path):
        masked_voronoi = MaskVoronoiDiagram(voronoi, centerline_remaining)
        voronoi_remaining = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
        write_polydata(voronoi_remaining, voronoi_remaining_path)
    else:
        voronoi_remaining = read_polydata(voronoi_remaining_path)

    # Extract translation vectors
    print("Computing translation directions.")
    direction = "horizont"
    middle_points, middle_ids = get_spline_points(centerlines, beta, direction, region_points_vtk)
    dx_p1 = middle_points[0] - p1

    # Move centerline manually for updating inlet
    # and outlet positions used to compute
    # new centerlines
    if beta != 0.0:
        new_centerlines = move_centerlines(centerlines, dx_p1, p1, p2, diverging_id,
                                           diverging_centerlines, direction)

        # Move anterior bend horizontally.
        # Iterate over points P from Voronoi diagram and manioulate
        print("Adjusting Voronoi diagram")
        voronoi_remaining = move_voronoi_horizontally(dx_p1, voronoi_remaining,
                                                      centerline_remaining, id1, id2,
                                                      diverging_id, clip=False)
        voronoi_siphon = move_voronoi_horizontally(dx_p1, voronoi_siphon,
                                                   centerline_siphon, id1, id2,
                                                   diverging_id, clip=True,
                                                   diverging_centerline_ispresent=diverging_centerline_ispresent)
    else:
        if diverging_centerline_ispresent:
            new_centerlines = merge_data([centerlines, extract_single_line(diverging_centerlines, 0)])
        else:
            new_centerlines = centerlines

        new_surface = surface
        write_polydata(new_centerlines, new_centerlines_path_tmp)

    if alpha == 0.0 and beta != 0.0:
        new_voronoi = merge_data([voronoi_remaining, voronoi_siphon])
        new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    elif alpha != 0.0:
        # Vertical movement
        print("Moving geometry vertically")
        new_surface, new_centerlines = move_vessel_vertically(alpha, voronoi_remaining,
                                                              voronoi_siphon,
                                                              new_centerlines, region_points, poly_ball_size)

    print("Smoothing, clean, and check surface")
    write_polydata(new_centerlines, "test_new_cl_alpha_{}_beta_{}.vtp".format(alpha, beta))
    new_surface = prepare_surface_output(new_surface, surface,
                                         new_centerlines, output_filepath,
                                         test_merge=True, changed=True,
                                         old_centerline=merge_data([centerlines,
                                                                    diverging_centerlines]))
    write_polydata(new_centerlines, new_centerlines_path)
    write_polydata(new_surface, output_filepath)


def move_vessel_vertically(alpha, voronoi_remaining,
                           voronoi_siphon, centerlines, region_points, poly_ball_size):
    """
    Secondary script used for vertical displacement of
    the blood vessel. Moves the input voronoi diagram and
    centerline in the vertical direction.

    Args:
        centerlines (vtkPolyData): Centerline through the geometry.
        alpha (float): Extension / Compression factor in vertical direction.
        voronoi_remaining (vtkPolyData): Voronoi diagram excluding siphon.
        voronoi_siphon (vtkPolyData): Voronoi diagram representing siphon.
        region_points (list): Points defining the siphon to be manipulated.
        poly_ball_size (list): Resulution of surface model.

    Returns:
        new_surface (vtkPolyData): New surface model.
    Returns:
        new_centerline (vtkPolyData): New centerline.
    """

    # Set and get clipping points, centerlines and diverging centerlines
    centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        find_region_of_interest_and_diverging_centerlines(centerlines, region_points)
    diverging_centerline_ispresent = diverging_centerlines is not None
    diverging_id = None if len(diverging_ids) == 0 else diverging_ids[0]

    # Special cases including the ophthalmic artery
    if diverging_centerline_ispresent:
        patch_diverging_line = clip_diverging_line(extract_single_line(diverging_centerlines, 0), region_points[0],
                                                   diverging_id)

    # Get clipped curve
    print("Clipping centerlines.")
    locator = get_locator(extract_single_line(centerlines, 0))
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    p1 = centerlines.GetPoint(id1)
    p2 = centerlines.GetPoint(id2)
    centerline_siphon = extract_single_line(centerlines, 0, startID=id1, endID=id2)

    if diverging_centerline_ispresent:
        eyeline_end = extract_single_line(patch_diverging_line, 1)
        centerline_siphon = merge_data([centerline_siphon, eyeline_end])

    # Find ID of middle pooint:
    print("Finding points to spline through.")
    direction = "vertical"
    middle_points, middle_ids, dx = get_spline_points(extract_single_line(centerlines, 0), alpha, direction,
                                                      region_points_vtk)

    # Iterate over points P from Voronoi diagram and manipulate
    print("Adjust voronoi diagram")
    voronoi_siphon = move_voronoi_vertically(voronoi_siphon, centerline_siphon, id1,
                                             diverging_id, dx, diverging_centerline_ispresent)
    new_voronoi = merge_data([voronoi_remaining, voronoi_siphon])

    # Move centerline manually for postprocessing
    new_centerlines = move_centerlines(centerlines, dx, p1, p2, diverging_id, diverging_centerlines, direction)

    # Write a new surface from the new voronoi diagram
    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    return new_surface, new_centerlines


def move_voronoi_horizontally(dx_p1, voronoi_clipped, centerline_clipped, id1, id2,
                              clip_id, clip=False, diverging_centerline_ispresent=False):
    """
    Iterate through voronoi diagram and move based on a profile
    for horizontal movement. Includes special treatment of
    opthalmic artery if present.

    Args:
        dx_p1 (ndarray): Direction to move upstream.
        voronoi_clipped (vtkPolyData): Voronoi diagram to be moved.
        centerline_clipped (vtkPolyData): Centerline corresponding voronoi diagram.
        id1 (int): Index of first clipping point.
        id2 (int): Index of second clipping point.
        clip_id (int): Index where opthamlic artery is located (if present)
        clip (bool): Determines which part of geometry is being moved, True if siphon.
        diverging_centerline_ispresent (bool): Determines presence of opthamlic artery.

    Returns:
        newDataSet (vtkPolyData): Manipulated Voronoi diagram.
    """

    centerline_loc = get_locator(centerline_clipped)
    newDataSet = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if clip:
        # Find boundaries
        idmid_0 = int((id1 + id2) / 2.)
        id1_0 = id1
        id1 = 0
        if diverging_centerline_ispresent:
            l1 = extract_single_line(centerline_clipped, 0)
            id2 = len(get_curvilinear_coordinate(l1))
        else:
            id2 = len(get_curvilinear_coordinate(centerline_clipped))
        idmid = int((id1 + id2) / 2.)

        # Manpipulation of voronoi diagram..
        if diverging_centerline_ispresent:
            # ..with opthalmic artery
            for p in range(voronoi_clipped.GetNumberOfPoints()):
                cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

                if cl_id < idmid:
                    dist = dx_p1 * (idmid ** 2 - cl_id ** 2) / (idmid ** 2 - id1 ** 2)
                else:
                    if cl_id <= (id2 - 1):
                        dist = -dx_p1 * (cl_id - idmid) ** 0.5 / (id2 - idmid) ** 0.5
                    else:
                        # Move opthalmic artery based on its discovery ID
                        if clip_id < idmid_0:
                            cl_id = clip_id - id1_0
                            dist = dx_p1 * (idmid ** 2 - cl_id ** 2) / (idmid ** 2 - id1 ** 2)
                        else:
                            cl_id = clip_id - id1_0
                            dist = -dx_p1 * (cl_id - idmid) ** 0.5 / (id2 - idmid) ** 0.5

                points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
                verts.InsertNextCell(1)
                verts.InsertCellPoint(p)

        else:
            # ..without opthalmic artery
            for p in range(voronoi_clipped.GetNumberOfPoints()):
                cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

                if cl_id < idmid:
                    dist = dx_p1 * (idmid ** 2 - cl_id ** 2) / (idmid ** 2 - id1 ** 2)
                else:
                    dist = -dx_p1 * (cl_id - idmid) ** 0.5 / (id2 - idmid) ** 0.5

                points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
                verts.InsertNextCell(1)
                verts.InsertCellPoint(p)

    else:
        # Move reamining part of the voronoi diagram
        # representing the geometry excluding the siphon to be moved
        for p in range(voronoi_clipped.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

            if cl_id <= id1:
                dist = dx_p1
            else:
                dist = -dx_p1

            points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    # Insert points and create new voronoi diagram
    newDataSet.SetPoints(points)
    newDataSet.SetVerts(verts)
    newDataSet.GetPointData().AddArray(voronoi_clipped.GetPointData().GetArray(radiusArrayName))

    return newDataSet


def move_voronoi_vertically(voronoi_clipped, centerline_clipped, id1_0, clip_id,
                            dx, diverging_centerline_ispresent=False):
    """
    Iterate through voronoi diagram and move based on a profile
    for vertical movement. Includes special treatment of
    opthalmic artery if present.

    Args:
        voronoi_clipped (vtkPolyData): Voronoi diagram to be moved.
        centerline_clipped (vtkPolyData): Centerline corresponding voronoi diagram.
        id1_0 (int): Index of first clipping point.
        clip_id (int): Index where opthamlic artery is located (if present)
        dx (ndarray): Direction to move.
        diverging_centerline_ispresent (bool): Determines presence of opthamlic artery.

    Returns:
        newDataSet (vtkPolyData): Manipulated Voronoi diagram.
    """

    centerline_loc = get_locator(centerline_clipped)
    newDataSet = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    # Manpipulation of voronoi diagram..
    if diverging_centerline_ispresent:
        # ..with opthalmic artery
        id1 = I1 = 0
        l1 = extract_single_line(centerline_clipped, 0)
        id2 = len(get_curvilinear_coordinate(l1))
        I2 = id2 - 1
        for p in range(voronoi_clipped.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

            if cl_id <= I2:
                dist = 4 * dx * (cl_id - I1) * (I2 - cl_id) / (I2 - I1) ** 2
            else:
                cl_id = clip_id - id1_0
                dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2

            points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    else:
        # ..witout opthalmic artery
        id1 = 0
        id2 = len(get_curvilinear_coordinate(centerline_clipped)) - 1
        for p in range(voronoi_clipped.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

            dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2

            points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    # Insert points and create new voronoi diagram
    newDataSet.SetPoints(points)
    newDataSet.SetVerts(verts)
    newDataSet.GetPointData().AddArray(voronoi_clipped.GetPointData().GetArray(radiusArrayName))

    return newDataSet


def read_command_line():
    """
    Read arguments from commandline
    """
    # Description of the script
    description = "Moves a selected part of a tubular geometry, " + \
                  "in two (horizontal, vertical) geometry-spesific directions. " + \
                  "Magnitude of movement is defined by the parameters alpha and beta" + \
                  "Primary script used for application in blood vessels."

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required named arguments')

    # Required arguments
    required.add_argument('-i', '--ifile', type=str, default=None,
                          help="Path to the surface model", required=True)
    required.add_argument("-o", "--ofile", type=str, default=None, required=True,
                          help="Relative path to the output surface. The default folder is" +
                               " the same as the input file, and a name with a combination of the" +
                               " parameters.")
    # Optional arguments
    parser.add_argument('-s', '--smooth', type=bool, default=True,
                        help="Smooth the voronoi diagram, default is True")
    parser.add_argument('-f', '--smooth_factor', type=float, default=0.25,
                        help="If smooth option is true then each voronoi point" +
                             " that has a radius less then MISR*(1-smooth_factor) at" +
                             " the closest centerline point is removed.")
    parser.add_argument("-n", "--no_smooth", type=bool, default=False,
                        help="If true and smooth is true the user, if no_smooth_point is" +
                             " not given, the user can provide points where the surface not will" +
                             " be smoothed.")
    parser.add_argument("--no_smooth_point", nargs="+", type=float, default=None,
                        help="If model is smoothed the user can manually select points on" +
                             " the surface that will not be smoothed. A centerline will be" +
                             " created to the extra point, and the section were the centerline" +
                             " differ from the other centerlines will be keept un-smoothed. This" +
                             " can be practicle for instance when manipulating geometries" +
                             " with aneurysms")
    parser.add_argument("-b", "--poly-ball-size", nargs=3, type=int, default=[120, 120, 120],
                        help="The size of the poly balls that will envelope the new" +
                             " surface. The default value is 120, 120, 120. If two tubular" +
                             " structures are very close compared to the bounds, the poly ball" +
                             " size should be adjusted", metavar="size")

    # Set region of interest:
    parser.add_argument("-r", "--region-of-interest", type=str, default="manuall",
                        choices=["manuall", "commandline", "landmarking"],
                        help="The method for defining the region to be changed. There are" +
                             " three options: 'manuall', 'commandline', 'landmarking'. In" +
                             " 'manuall' the user will be provided with a visualization of the" +
                             " input surface, and asked to provide an end and start point of the" +
                             " region of interest. Note that not all algorithms are robust over" +
                             " bifurcations. If 'commandline' is provided, then '--region-points'" +
                             " is expected to be provided. Finally, if 'landmarking' is" +
                             " given, it will look for the output from running" +
                             " automated_geometric_quantities.py.")
    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" +
                             " argument have to be given. The method expects two points" +
                             " which defines the start and end of the region of interest. If" +
                             " 'method' is set to stenosis, then one point can be provided as well," +
                             " which is assumbed to be the center of a new stenosis." +
                             " Example providing the points (1, 5, -1) and (2, -4, 3):" +
                             " --stenosis-points 1 5 -1 2 -4 3")
    # "Variation" argments
    parser.add_argument("--alpha", type=float, default=0.0,
                        help="Compression factor in vertical direction, " +
                             "ranging from -1.0 to 1.0, defining the magnitude " +
                             "of streching or compression of the tubular structure.")
    parser.add_argument("--beta", type=float, default=0.0,
                        help="Compression factor in verti cal direction,  " +
                             "ranging from -1.0 to 1.0, defining the magnitude " +
                             "of streching or compression of the tubular structure.")
    # Outputfile argument
    args = parser.parse_args()

    if args.no_smooth_point is not None and len(args.no_smooth_point):
        if len(args.no_smooth_point) % 3 != 0:
            raise ValueError("ERROR: Please provide the no smooth point(s) as a multiple of 3")

    return dict(input_filepath=args.ifile, smooth=args.smooth,
                smooth_factor=args.smooth_factor, alpha=args.alpha, beta=args.beta,
                output_filepath=args.ofile, poly_ball_size=args.poly_ball_size,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point,
                region_of_interest=args.region_of_interest, region_points=args.region_points)


if __name__ == "__main__":
    move_vessel(**read_command_line())

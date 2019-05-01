##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from morphman.common.argparse_common import *
from morphman.common.surface_operations import *


def manipulate_bend(input_filepath, output_filepath, smooth, smooth_factor, region_of_interest, region_points,
                    alpha, beta, poly_ball_size, no_smooth, no_smooth_point, resampling_step):
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

    Continuation in the manipulate_bend_vertically-method for vertical movement.

    Args:
        input_filepath (str): Path to input surface.
        output_filepath (str): Path to output the manipulated surface.
        smooth (bool): Smooth the Voronoi diagram.
        smooth_factor (float): Smoothing factor used for Voronoi diagram smoothing.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'landmarking']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        alpha (float): Extension / Compression factor in vertical direction.
        beta (float): Extension / Compression factor in horizontal direction.
        poly_ball_size (list): Resolution of polyballs used to create surface.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
        resampling_step (float): Resampling length when resampling centerline.
    """

    # Input filenames
    base_path = get_path_names(input_filepath)

    # Centerlines filenames
    centerline_clipped_path = base_path + "_centerline_clipped.vtp"
    centerline_clipped_part_path = base_path + "_centerline_clipped_part.vtp"
    new_centerlines_path = base_path + "_centerlines_alpha_%s_beta_%s.vtp" % (alpha, beta)
    new_centerlines_path_tmp = base_path + "_new_centerlines_alpha_%s_beta_%s_tmp.vtp" % (alpha, beta)
    centerlines_path = base_path + "_centerline.vtp"

    # Voronoi diagrams filenames
    voronoi_remaining_path = base_path + "_voronoi_bend_remaining.vtp"
    voronoi_bend_path = base_path + "_voronoi_bend.vtp"

    # Region point filename
    point_path = base_path + "_anterior_bend.particles"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Compute centerlines
    inlet, outlets = get_inlet_and_outlet_centers(surface, base_path)

    print("-- Compute centerlines and Voronoi diagram")
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=resampling_step,
                                                         smooth=False, base_path=base_path)
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids, resampling_step)

    # Get region of interest
    if region_of_interest == "landmarking":
        if not path.exists(point_path):
            raise RuntimeError(("The given .particles file: %s does not exist. Please run" +
                                " landmarking with automated_landmarking.py first.") % point_path)
        region_points = np.loadtxt(point_path)
    else:
        _, _, _, region_points, _ = get_line_to_change(capped_surface, centerlines,
                                                       region_of_interest, "bend", region_points, 0)
        region_points = [[region_points[3 * i], region_points[3 * i + 1], region_points[3 * i + 2]]
                         for i in range(len(region_points) // 3)]

    # Set and get clipping points, centerlines and diverging centerlines
    centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        get_region_of_interest_and_diverging_centerlines(centerlines, region_points)
    diverging_centerline_ispresent = diverging_centerlines is not None
    diverging_id = None if len(diverging_ids) == 0 else diverging_ids[0]

    # Handle diverging centerlines within region of interest
    if diverging_centerline_ispresent:
        print("-- Clipping a centerline divering in the region of interest.")
        patch_diverging_line = get_clipped_diverging_centerline(extract_single_line(diverging_centerlines, 0),
                                                                region_points[0], diverging_id)

    # Clip centerline
    print("-- Clipping centerlines.")
    locator = get_vtk_point_locator(extract_single_line(centerlines, 0))
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    p1 = centerlines.GetPoint(id1)
    p2 = centerlines.GetPoint(id2)
    centerline_remaining = create_parent_artery_patches(centerlines,
                                                        region_points_vtk, siphon=True)
    centerline_bend = extract_single_line(centerlines, 0, start_id=id1, end_id=id2)

    if diverging_centerline_ispresent:
        diverging_centerline_end = extract_single_line(patch_diverging_line, 1)
        centerline_bend = vtk_merge_polydata([centerline_bend, diverging_centerline_end])

    write_polydata(centerline_remaining, centerline_clipped_path)
    write_polydata(centerline_bend, centerline_clipped_part_path)

    # Clip Voronoi diagram into
    # bend and remaining part of geometry
    print("-- Clipping Voronoi diagrams")
    voronoi_bend, voronoi_remaining = get_split_voronoi_diagram(voronoi, [centerline_bend, centerline_remaining])
    write_polydata(voronoi_bend, voronoi_bend_path)
    write_polydata(voronoi_remaining, voronoi_remaining_path)

    # Extract translation vectors
    print("-- Computing translation directions.")
    direction = "horizont"
    middle_points, middle_ids = get_direction_parameters(centerlines, beta, direction, region_points_vtk)
    dx_p1 = middle_points[0] - p1

    # Move centerline manually for updating inlet
    # and outlet positions used to compute
    # new centerlines
    if beta != 0.0:
        new_centerlines = get_manipulated_centerlines(centerlines, dx_p1, p1, p2, diverging_id,
                                                      diverging_centerlines, direction)

        # Move anterior bend horizontally.
        # Iterate over points P from Voronoi diagram and manipulate
        print("-- Adjusting Voronoi diagram")
        voronoi_remaining = move_voronoi_horizontally(dx_p1, voronoi_remaining,
                                                      centerline_remaining, id1, id2,
                                                      diverging_id, clip=False)
        voronoi_bend = move_voronoi_horizontally(dx_p1, voronoi_bend,
                                                 centerline_bend, id1, id2,
                                                 diverging_id, clip=True,
                                                 diverging_centerline_ispresent=diverging_centerline_ispresent)
    else:
        if diverging_centerline_ispresent:
            new_centerlines = vtk_merge_polydata([centerlines, extract_single_line(diverging_centerlines, 0)])
        else:
            new_centerlines = centerlines

        new_surface = surface
        write_polydata(new_centerlines, new_centerlines_path_tmp)

    if alpha == 0.0 and beta != 0.0:
        print("-- Creating new surface.")
        new_voronoi = vtk_merge_polydata([voronoi_remaining, voronoi_bend])
        new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    elif alpha != 0.0:
        # Update the region points
        locator = get_vtk_point_locator(new_centerlines)
        region_points[0] = new_centerlines.GetPoint(locator.FindClosestPoint(region_points[0]))
        region_points[1] = new_centerlines.GetPoint(locator.FindClosestPoint(region_points[1]))

        # Vertical movement
        print("-- Moving geometry vertically")
        new_surface, new_centerlines = manipulate_bend_vertically(alpha, voronoi_remaining,
                                                                  voronoi_bend,
                                                                  new_centerlines, region_points, poly_ball_size)

    print("-- Smoothing, clean, and check surface")
    new_surface = prepare_output_surface(new_surface, surface,
                                         new_centerlines, output_filepath,
                                         test_merge=True, changed=True,
                                         old_centerline=vtk_merge_polydata([centerlines,
                                                                            diverging_centerlines]))
    write_polydata(new_centerlines, new_centerlines_path)
    write_polydata(new_surface, output_filepath)


def manipulate_bend_vertically(alpha, voronoi_remaining, voronoi_bend, centerlines, region_points, poly_ball_size):
    """
    Secondary script used for vertical displacement of
    the blood vessel. Moves the input voronoi diagram and
    centerline in the vertical direction.

    Args:
        centerlines (vtkPolyData): Centerline through the geometry.
        alpha (float): Extension / Compression factor in vertical direction.
        voronoi_remaining (vtkPolyData): Voronoi diagram excluding bend.
        voronoi_bend (vtkPolyData): Voronoi diagram representing bend.
        region_points (list): Points defining the bend to be manipulated.
        poly_ball_size (list): Resolution of surface model.

    Returns:
        new_surface (vtkPolyData): New surface model.
    Returns:
        new_centerline (vtkPolyData): New centerline.
    """

    # Set and get clipping points, centerlines and diverging centerlines
    centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        get_region_of_interest_and_diverging_centerlines(centerlines, region_points)
    diverging_centerline_ispresent = diverging_centerlines is not None
    diverging_id = None if len(diverging_ids) == 0 else diverging_ids[0]

    # Special cases including the diverging centerline
    if diverging_centerline_ispresent:
        patch_diverging_line = get_clipped_diverging_centerline(extract_single_line(diverging_centerlines, 0),
                                                                region_points[0], diverging_id)

    # Get clipped curve
    print("-- Clipping centerlines.")
    locator = get_vtk_point_locator(extract_single_line(centerlines, 0))
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    p1 = centerlines.GetPoint(id1)
    p2 = centerlines.GetPoint(id2)
    centerline_bend = extract_single_line(centerlines, 0, start_id=id1, end_id=id2)

    if diverging_centerline_ispresent:
        diverging_centerline_end = extract_single_line(patch_diverging_line, 1)
        centerline_bend = vtk_merge_polydata([centerline_bend, diverging_centerline_end])

    # Find ID of middle point:
    print("-- Finding horizontal direction parameters.")
    direction = "vertical"
    middle_points, middle_ids, dx = get_direction_parameters(extract_single_line(centerlines, 0), alpha, direction,
                                                             region_points_vtk)

    # Iterate over points P from Voronoi diagram and manipulate
    print("-- Adjust Voronoi diagram")
    voronoi_bend = move_voronoi_vertically(voronoi_bend, centerline_bend, id1,
                                           diverging_id, dx, diverging_centerline_ispresent)
    new_voronoi = vtk_merge_polydata([voronoi_remaining, voronoi_bend])

    # Move centerline manually for postprocessing
    new_centerlines = get_manipulated_centerlines(centerlines, dx, p1, p2, diverging_id, diverging_centerlines,
                                                  direction)

    # Write a new surface from the new voronoi diagram
    print("-- Creating new surface.")
    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    return new_surface, new_centerlines


def move_voronoi_horizontally(dx_p1, voronoi_clipped, centerline_clipped, id1, id2,
                              clip_id, clip=False, diverging_centerline_ispresent=False):
    """
    Iterate through voronoi diagram and move based on a profile
    for horizontal movement. Includes special treatment of
    diverging centerlines if present.

    Args:
        dx_p1 (ndarray): Direction to move upstream.
        voronoi_clipped (vtkPolyData): Voronoi diagram to be moved.
        centerline_clipped (vtkPolyData): Centerline corresponding voronoi diagram.
        id1 (int): Index of first clipping point.
        id2 (int): Index of second clipping point.
        clip_id (int): Index where diverging centerline is located (if present)
        clip (bool): Determines which part of geometry is being moved, True if bend.
        diverging_centerline_ispresent (bool): Determines presence of diverging centerline.

    Returns:
        new_dataset (vtkPolyData): Manipulated Voronoi diagram.
    """

    centerline_loc = get_vtk_point_locator(centerline_clipped)
    new_dataset = vtk.vtkPolyData()
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

        # Manipulation of voronoi diagram..
        if diverging_centerline_ispresent:
            # ..with diverging centerline
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
            # ..without diverging centerline
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
        # Move remaining part of the voronoi diagram
        # representing the geometry excluding the bend to be moved
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
    new_dataset.SetPoints(points)
    new_dataset.SetVerts(verts)
    new_dataset.GetPointData().AddArray(voronoi_clipped.GetPointData().GetArray(radiusArrayName))

    return new_dataset


def move_voronoi_vertically(voronoi_clipped, centerline_clipped, id1_0, clip_id,
                            dx, diverging_centerline_ispresent=False):
    """
    Iterate through voronoi diagram and move based on a profile
    for vertical movement. Includes special treatment of
    diverging centerline if present.

    Args:
        voronoi_clipped (vtkPolyData): Voronoi diagram to be moved.
        centerline_clipped (vtkPolyData): Centerline corresponding voronoi diagram.
        id1_0 (int): Index of first clipping point.
        clip_id (int): Index where diverging centerline is located (if present)
        dx (ndarray): Direction to move.
        diverging_centerline_ispresent (bool): Determines presence of diverging centerline.

    Returns:
        new_dataset (vtkPolyData): Manipulated Voronoi diagram.
    """

    centerline_loc = get_vtk_point_locator(centerline_clipped)
    new_dataset = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    # Manipulation of voronoi diagram..
    if diverging_centerline_ispresent:
        # ..with diverging centerline
        l1 = extract_single_line(centerline_clipped, 0)
        id1 = 0
        id2 = len(get_curvilinear_coordinate(l1)) - 1
        for p in range(voronoi_clipped.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

            if cl_id <= id2:
                dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2
            else:
                cl_id = clip_id - id1_0
                dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2

            points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    else:
        # ..without diverging centerline
        id1 = 0
        id2 = len(get_curvilinear_coordinate(centerline_clipped)) - 1
        for p in range(voronoi_clipped.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(voronoi_clipped.GetPoint(p))

            dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2

            points.InsertNextPoint(np.asarray(voronoi_clipped.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    # Insert points and create new voronoi diagram
    new_dataset.SetPoints(points)
    new_dataset.SetVerts(verts)
    new_dataset.GetPointData().AddArray(voronoi_clipped.GetPointData().GetArray(radiusArrayName))

    return new_dataset


def read_command_line_bend(input_path=None, output_path=None):
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
    parser.add_argument("--alpha", type=float, default=0.0,
                        help="Compression factor in vertical direction, " +
                             "ranging from -1.0 to 1.0, defining the magnitude " +
                             "of stretching or compression of the tubular structure.")
    parser.add_argument("--beta", type=float, default=0.0,
                        help="Compression factor in vertical direction,  " +
                             "ranging from -1.0 to 1.0, defining the magnitude " +
                             "of stretching or compression of the tubular structure.")
    # Output file argument
    if required:
        args = parser.parse_args()
    else:
        args = parser.parse_args(["-i" + input_path, "-o" + output_path])

    if args.no_smooth_point is not None and len(args.no_smooth_point):
        if len(args.no_smooth_point) % 3 != 0:
            raise ValueError("ERROR: Please provide the no smooth point(s) as a multiple of 3")

    return dict(input_filepath=args.ifile, smooth=args.smooth,
                smooth_factor=args.smooth_factor, alpha=args.alpha, beta=args.beta,
                output_filepath=args.ofile, poly_ball_size=args.poly_ball_size,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point,
                region_of_interest=args.region_of_interest, region_points=args.region_points,
                resampling_step=args.resampling_step)


def main_bend():
    manipulate_bend(**read_command_line_bend())


if __name__ == "__main__":
    manipulate_bend(**read_command_line_bend())

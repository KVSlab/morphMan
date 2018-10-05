##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from common import *
from argparse_common import *


def curvature_variations(input_filepath, smooth, smooth_factor, smooth_factor_line, iterations,
                         smooth_line, output_filepath, poly_ball_size, region_of_interest,
                         region_points, resampling_step, no_smooth, no_smooth_point):
    """
    Create a sharper or smoother version of the input geometry,
    determined by a smoothed version of the siphon centerline.

    Args:
        input_filepath (str): Path to input surface.
        output_filepath (str): Path to output the manipulated surface.
        smooth (bool): Smooth the Voronoi diagram.
        smooth_line (bool): Smooth centerline if True, anti-smooth if False.
        smooth_factor (float): Smoothing factor used for Voronoi diagram smoothing.
        smooth_factor_line (float): Smoothing factor used for centerline smoothing.
        iterations (int): Smoothing iterations of centerline.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'first_line']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        poly_ball_size (list): Resolution of polyballs used to create surface.
        resampling_step (float): Resampling length for centerline resampling.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
    """
    # Input filenames
    base_path = get_path_names(input_filepath)

    # Centerlines
    centerlines_path = base_path + "_centerline.vtp"
    centerline_smooth_path = base_path + "_centerline_smoothed.vtp"
    new_centerlines_path = base_path + "_centerline_new.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Voronoi diagrams filenames
    voronoi_remaining_path = base_path + "_voronoi_remaining.vtp"
    voronoi_region_path = base_path + "_voronoi_region.vtp"
    voronoi_diverging_path = base_path + "_voronoi_diverging.vtp"

    # Compute centerlines
    inlet, outlets = get_centers(surface, base_path)

    print("-- Compute centerlines and Voronoi diagram")
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=resampling_step,
                                                         smooth=False, base_path=base_path)
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)
    # Get region of interest
    _, region_points = get_line_to_change(capped_surface, centerlines,
                                          region_of_interest, "variation", region_points, 0)
    region_points = [[region_points[3 * i], region_points[3 * i + 1], region_points[3 * i + 2]]
                     for i in range(len(region_points) // 3)]

    # Set and get clipping points, centerlines and diverging centerlines
    centerlines_complete, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        find_region_of_interest_and_diverging_centerlines(centerlines, region_points)
    diverging_centerline_ispresent = diverging_centerlines is not None
    diverging_id = None if len(diverging_ids) == 0 else diverging_ids[0]

    # Handle diverging centerlines within region of interest
    diverging_centerlines_patch = []
    if diverging_centerline_ispresent:
        print("-- Clipping a centerline divering in the region of interest.")
        for i in range(diverging_centerlines.GetNumberOfCells()):
            diverging_centerline = extract_single_line(diverging_centerlines, i)
            patch_diverging_centerline = clip_diverging_line(diverging_centerline,
                                                             region_points[0], diverging_id)
            diverging_centerlines_patch.append(extract_single_line(patch_diverging_centerline, 1))
        diverging_centerlines_patch = merge_data(diverging_centerlines_patch)

    # Clip centerline
    print("-- Clipping centerlines.")
    locator = get_locator(extract_single_line(centerlines_complete, 0))
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    centerline_remaining = create_parent_artery_patches(centerlines_complete,
                                                        region_points_vtk, siphon=True)
    centerline_region = extract_single_line(centerlines_complete, 0, startID=id1, endID=id2)

    if diverging_centerline_ispresent:
        centerline_region = merge_data([centerline_region, diverging_centerlines_patch])

    # Clip Voronoi diagram into
    # bend and remaining part of geometry
    print("-- Clipping Voronoi diagrams")
    voronoi_region, voronoi_remaining = split_voronoi_with_centerlines(voronoi,
                                                                       centerline_region,
                                                                       centerline_remaining)
    # Spearate diverging parts and main region
    if diverging_centerline_ispresent:
        voronoi_region, voronoi_diverging = split_voronoi_with_centerlines(voronoi_region,
                                                                           extract_single_line(centerlines_complete,
                                                                                               0,
                                                                                               startID=id1,
                                                                                               endID=id2),
                                                                           diverging_centerlines_patch)
        write_polydata(voronoi_diverging, voronoi_diverging_path)

    write_polydata(voronoi_region, voronoi_region_path)
    write_polydata(voronoi_remaining, voronoi_remaining_path)

    print("-- Smooth / sharpen centerline")
    smoothed_centerline_region = vmtk_centerline_geometry(centerline_region, True, True,
                                                          factor=smooth_factor_line, iterations=iterations)
    write_polydata(smoothed_centerline_region, centerline_smooth_path)

    print("-- Smooth / sharpen Voronoi diagram")
    moved_voronoi_region = make_voronoi_smooth(voronoi_region, centerline_region, smoothed_centerline_region,
                                               smooth_line)
    # Move diverging centerlines and combine all diagrams
    smoothed_centerline_region_siphon = extract_single_line(smoothed_centerline_region, 0)
    if diverging_centerline_ispresent:
        centerline_region_siphon = extract_single_line(centerline_region, 0)
        moved_voronoi_diverging = make_voronoi_smooth(voronoi_diverging, centerline_region_siphon,
                                                      smoothed_centerline_region_siphon, smooth_line,
                                                      div=True, div_point=diverging_centerlines.GetPoint(diverging_id))
        new_voronoi = merge_data([moved_voronoi_region, moved_voronoi_diverging, voronoi_remaining])
    else:
        new_voronoi = merge_data([moved_voronoi_region, voronoi_remaining])

    print("-- Moving centerlines")
    new_centerlines = move_all_centerlines(centerlines_complete, smoothed_centerline_region_siphon, diverging_id,
                                           diverging_centerlines, smooth_line)

    # Create new surface and move centerlines (for postprocessing)
    print("-- Create new surface")
    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    print("-- Smoothing, clean, and check surface")
    new_surface = prepare_surface_output(new_surface, surface,
                                         new_centerlines, output_filepath,
                                         test_merge=True, changed=True,
                                         old_centerline=centerlines_complete)

    write_polydata(new_centerlines, new_centerlines_path)
    write_polydata(new_surface, output_filepath)


def make_voronoi_smooth(voronoi, old_cl, new_cl, smooth_line, div=False, div_point=None):
    """
    Move the voronoi diagram based on a smoothed
    version of the centerline, returning a
    manipulated version of the voronoi diagram.

    Args:
        voronoi (vtkPolyData): Voronoi diagram data set.
        old_cl (vtkPolyData): Unsmoothed centerline points.
        new_cl (vtkPolyData): Smoothed centerline points.
        smooth_line (bool): Determines if model becomes smoother or sharper.
        div (bool): True if centerline is a diverging line.
        div_point (ndarray): Diverging point along siphon.

    Returns:
        new_dataset (vtkPolyData): Manipulated voronoi diagram.
    """
    locator = get_locator(old_cl)
    n = voronoi.GetNumberOfPoints()
    new_dataset = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    # Define segments for transitioning
    id_end = old_cl.GetNumberOfPoints() - 1
    id_midend = int(id_end * 0.9)
    id_startmid = int(id_end * 0.1)
    id_mid = int(id_end * 0.2)
    # Iterate through voronoi points
    for i in range(n):
        if div:
            cl_id = locator.FindClosestPoint(div_point)
        else:
            cl_id = locator.FindClosestPoint(voronoi.GetPoint(i))
        p0 = np.asarray(old_cl.GetPoint(cl_id))
        p1 = np.asarray(new_cl.GetPoint(cl_id))
        if smooth_line:
            dx = p1 - p0
        else:
            dx = -(p1 - p0)

        # Smooth transition at inlet and at end of siphon
        if not div:
            if cl_id < id_startmid:
                dx = 0
            elif id_startmid <= cl_id < id_mid:
                dx = dx * (cl_id - id_startmid) / float(id_mid - id_startmid)
            elif cl_id > id_midend:
                dx = dx * (id_end - cl_id) / float(id_end - id_midend)
        points.InsertNextPoint(np.asarray(voronoi.GetPoint(i)) + dx)
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)

    new_dataset.SetPoints(points)
    new_dataset.SetVerts(verts)
    new_dataset.GetPointData().AddArray(voronoi.GetPointData().GetArray(radiusArrayName))

    return new_dataset


def move_all_centerlines(old_cl, new_cl, diverging_id, diverging_centerlines, smooth_line):
    """
    Takes the original centerline and eventual diverging centerlines, performing
    manual manipulation by smoothing / sharpening the centerline based on a
    smoothed region of the original centerline.
    Returns the final complete and manipulated centerline.

    Args:
        old_cl (vtkPolyData): Centerlines excluding diverging centerlines.
        new_cl (vtkPolyData): Smoothed region of the centerline.
        diverging_id (int): List of index where centerlines diverge from region of interest.
        diverging_centerlines (vtkPolyData): Centerlines which diverge from region of interest.
        smooth_line (bool): Determines if model becomes smoother or sharper.

    Returns:
        centerline (vtkPolyData): Manipulated centerline.
    """
    if diverging_id is not None:
        old_cl = merge_data([old_cl, diverging_centerlines])

    number_of_points = old_cl.GetNumberOfPoints()
    number_of_cells = old_cl.GetNumberOfCells()

    centerline = vtk.vtkPolyData()
    centerline_points = vtk.vtkPoints()
    centerline_cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, number_of_points)

    count = 0
    id_end = new_cl.GetNumberOfPoints() - 1
    p1 = new_cl.GetPoint(0)
    p2 = new_cl.GetPoint(id_end)

    for i in range(number_of_cells):
        line = extract_single_line(old_cl, i)
        locator = get_locator(line)
        id1 = locator.FindClosestPoint(p1)
        if i == (number_of_cells - 1) and diverging_id is not None:
            id2 = diverging_id
        else:
            id2 = locator.FindClosestPoint(p2)
        n = line.GetNumberOfPoints()
        centerline_cell_array.InsertNextCell(n)
        radius_array_data = line.GetPointData().GetArray(radiusArrayName).GetTuple1

        # Iterate through voronoi points
        k = 0
        for p in range(n):
            p_old = np.asarray(line.GetPoint(p))
            if id1 < p < id2:
                p_new = np.asarray(new_cl.GetPoint(k))
                k += 1
                dx = (p_new - p_old)
            elif i == (number_of_cells - 1) and diverging_id is not None and p >= diverging_id:
                p_div = np.asarray(line.GetPoint(diverging_id))
                p_new = np.asarray(new_cl.GetPoint(k))
                dx = (p_new - p_div)
            else:
                dx = 0

            if not smooth_line:
                dx = - dx

            centerline_points.InsertNextPoint(p_old + dx)
            radius_array.SetTuple1(count, radius_array_data(p))
            centerline_cell_array.InsertCellPoint(count)
            count += 1

    centerline.SetPoints(centerline_points)
    centerline.SetLines(centerline_cell_array)
    centerline.GetPointData().AddArray(radius_array)

    return centerline


def read_command_line():
    """
    Read arguments from commandline
    """
    # Description of the script
    description = "Manipulates a selected part of a tubular geometry" + \
                  "by creaing a sharper or smoother version of the input geometry, " + \
                  "determined by a smoothed version of the centerline representing the selected part."

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required named arguments')

    # Add common arguments
    add_common_arguments(parser)

    # Set region of interest:
    parser.add_argument("-r", "--region-of-interest", type=str, default="manual",
                        choices=["manual", "commandline", "first_line"],
                        help="The method for defining the region to be changed. There are" +
                             " three options: 'manual', 'commandline', 'first_line'. In" +
                             " 'manual' the user will be provided with a visualization of the" +
                             " input surface, and asked to provide an end and start point of the" +
                             " region of interest. Note that not all algorithms are robust over" +
                             " bifurcations. If 'commandline' is provided, then '--region-points'" +
                             " is expected to be provided. Finally, if 'first_line' is" +
                             " given, the method will chose the first line of the geometry.")
    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" +
                             " argument have to be given. The method expects two points" +
                             " which defines the start and end of the region of interest. ")

    # Set problem spesific arguments
    parser.add_argument("-fl", "--smooth_factor_line", type=float, default=1.0,
                        help="Smoothing factor of centerline curve.")
    parser.add_argument("-it", "--iterations", type=int, default=100,
                        help="Smoothing iterations of centerline curve.")
    parser.add_argument("-sl", "--smooth_line", type=bool, default=True,
                        help="Smoothes centerline if True, anti-smoothes if False")

    # Parse
    args = parser.parse_args()

    return dict(input_filepath=args.ifile, smooth=args.smooth,
                smooth_factor=args.smooth_factor, smooth_factor_line=args.smooth_factor_line,
                iterations=args.iterations, smooth_line=args.smooth_line, output_filepath=args.ofile,
                poly_ball_size=args.poly_ball_size, region_of_interest=args.region_of_interest,
                region_points=args.region_points, resampling_step=args.resampling_step,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point)


if __name__ == "__main__":
    curvature_variations(**read_command_line())

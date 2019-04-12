##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from morphman.common.argparse_common import *
from morphman.common.surface_operations import *


def manipulate_curvature(input_filepath, smooth, smooth_factor, smooth_factor_line, iterations,
                         smooth_line, output_filepath, poly_ball_size, region_of_interest,
                         region_points, resampling_step, no_smooth, no_smooth_point):
    """
    Create a sharper or smoother version of the input geometry,
    determined by a smoothed version of the centerline.

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
    direction = "smoothed" if smooth_line else "extended"
    centerlines_path = base_path + "_centerline.vtp"
    centerline_smooth_path = base_path + "_centerline_smoothed.vtp"
    centerline_spline_path = base_path + "_centerline_splined.vtp"
    centerline_diverging_path = base_path + "_centerline_diverging.vtp"
    centerline_remaining_path = base_path + "_centerline_remaining.vtp"
    new_centerlines_path = base_path + "_centerline_new_%s.vtp" % direction

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Voronoi diagrams filenames
    voronoi_remaining_path = base_path + "_voronoi_curvature_remaining.vtp"
    voronoi_region_path = base_path + "_voronoi_region.vtp"
    voronoi_div_path = base_path + "_voronoi_diverging_{}.vtp"

    # Compute centerlines
    inlet, outlets = get_inlet_and_outlet_centers(surface, base_path)
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=resampling_step,
                                                         smooth=False, base_path=base_path)
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)
    # Get region of interest
    centerline_splined, centerline_remaining, \
            centerline_diverging, region_points, diverging_ids = get_line_to_change(capped_surface, centerlines,
                                                                                    region_of_interest, "variation",
                                                                                    region_points,
                                                                                    None)
    write_polydata(centerline_splined, centerline_spline_path)
    write_polydata(centerline_remaining, centerline_remaining_path)
    if centerline_diverging is not None:
        write_polydata(merge_data(centerline_diverging), centerline_diverging_path)

    # Split the Voronoi diagram
    print("-- Clipping Voronoi diagram")
    centerline_regions = [centerline_splined, centerline_remaining]
    if centerline_diverging is not None:
        for i, div_cl in enumerate(centerline_diverging):
            centerline_regions += [extract_single_line(div_cl, 0,
                                                       startID=diverging_ids[i][1])]
    voronoi_regions = split_voronoi_with_centerlines(voronoi, centerline_regions)

    write_polydata(voronoi_regions[0], voronoi_region_path)
    write_polydata(voronoi_regions[1], voronoi_remaining_path)
    for i in range(2, len(voronoi_regions)):
        write_polydata(voronoi_regions[i], voronoi_div_path.format(i-1))

    # Move the centerline
    print("-- Smooth / sharpen centerline")
    smoothed_centerline_splined = vmtk_centerline_geometry(centerline_splined, True, True,
                                                           factor=smooth_factor_line, iterations=iterations)
    write_polydata(smoothed_centerline_splined, centerline_smooth_path)

    # Move the Voronoi diagram
    print("-- Smooth / sharpen Voronoi diagram")
    diverging_points = [centerline_diverging[i].GetPoint(diverging_ids[i][1]) for i in range(len(diverging_ids))]
    moved_voronoi_region = make_voronoi_smooth(voronoi_regions[0], centerline_splined,
                                               smoothed_centerline_splined,
                                               smooth_line, voronoi_regions[2:],
                                               diverging_points)
    new_voronoi = merge_data([voronoi_regions[1]] + moved_voronoi_region)

    print("-- Moving centerlines")
    new_centerlines = move_all_centerlines(centerlines, smoothed_centerline_splined,
                                           diverging_ids, centerline_diverging, smooth_line)
    write_polydata(new_centerlines, new_centerlines_path)

    # Create new surface and move centerlines (for postprocessing)
    print("-- Create new surface")
    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    print("-- Smoothing, clean, and check surface")
    new_surface = prepare_output_surface(new_surface, surface,
                                         new_centerlines, output_filepath,
                                         test_merge=True, changed=True,
                                         old_centerline=centerlines)

    write_polydata(new_surface, output_filepath)


def make_voronoi_smooth(voronoi, old_cl, new_cl, smooth_line, div_voronoi, div_points):
    """
    Move the voronoi diagram based on a smoothed
    version of the centerline, returning a
    manipulated version of the voronoi diagram.

    Args:
        voronoi (vtkPolyData): Voronoi diagram data set.
        old_cl (vtkPolyData): Unsmoothed centerline points.
        new_cl (vtkPolyData): Smoothed centerline points.
        smooth_line (bool): Determines if model becomes smoother or sharper.
        div_voronoi (list): A list of diverging centerlines.
        div_points (list): List of diverging points along the region of interest.

    Returns:
        new_dataset (vtkPolyData): Manipulated voronoi diagram.
    """
    locator = vtk_point_locator(old_cl)
    n = voronoi.GetNumberOfPoints()
    new_dataset = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    # Define segments for transitioning
    id_end = old_cl.GetNumberOfPoints() - 1
    id_mid = int(id_end * 0.9)
    id_start = int(id_end * 0.1)

    # Iterate through voronoi points
    for i in range(n):
        cl_id = locator.FindClosestPoint(voronoi.GetPoint(i))
        p0 = np.asarray(old_cl.GetPoint(cl_id))
        p1 = np.asarray(new_cl.GetPoint(cl_id))
        if smooth_line:
            dx = p1 - p0
        else:
            dx = -(p1 - p0)

        # Smooth transition at inlet and at end of region of interest
        if cl_id < id_start:
            dx = dx * cl_id / float(id_start)
        elif cl_id > id_mid:
            dx = dx * (id_end - cl_id) / float(id_end - id_mid)

        points.InsertNextPoint(np.asarray(voronoi.GetPoint(i)) + dx)
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)

    new_dataset.SetPoints(points)
    new_dataset.SetVerts(verts)
    new_dataset.GetPointData().AddArray(voronoi.GetPointData().GetArray(radiusArrayName))

    # Offset diverging centerlines
    if div_voronoi != []:
        for i in range(len(div_voronoi)):
            cl_id = locator.FindClosestPoint(div_points[i])
            p0 = np.asarray(old_cl.GetPoint(cl_id))
            p1 = np.asarray(new_cl.GetPoint(cl_id))

            if smooth_line:
                dx = p1 - p0
            else:
                dx = -(p1 - p0)

            # Smooth transition at inlet and at end of region of interest
            if cl_id < id_start:
                dx = dx * cl_id / float(id_start)
            elif cl_id > id_mid:
                dx = dx * (id_end - cl_id) / float(id_end - id_mid)

            div_voronoi[i] = offset_voronoi(div_voronoi[i], dx)

    return [new_dataset] + div_voronoi


def move_all_centerlines(old_cl, new_cl, diverging_id, diverging_centerlines, smooth_line):
    """
    Takes the original centerline and eventual diverging centerlines, performing
    manual manipulation by smoothing / sharpening the centerline based on a
    smoothed region of the original centerline.
    Returns the final complete and manipulated centerline.

    Args:
        old_cl (vtkPolyData): Centerlines excluding diverging centerlines.
        new_cl (vtkPolyData): Smoothed region of the centerline.
        diverging_id (list): List of index where centerlines diverge from region of interest.
        diverging_centerlines (list): Centerlines which diverge from region of interest.
        smooth_line (bool): Determines if model becomes smoother or sharper.

    Returns:
        centerline (vtkPolyData): Manipulated centerline.
    """
    if diverging_id is not None:
        old_cl = vtk_merge_polydata([old_cl, diverging_centerlines])

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

    tol = get_tolerance(old_cl)

    # Iterate through centerline points
    div_count = -1
    for i in range(number_of_cells):
        line = extract_single_line(old_cl, i)

        # Check if line goes through the region of interest
        locator = vtk_point_locator(line)
        id1 = locator.FindClosestPoint(p1)
        id2 = locator.FindClosestPoint(p2)
        in_p1 = distance(line.GetPoint(id1), p1) < tol * 3
        in_p2 = distance(line.GetPoint(id2), p2) < tol * 3

        # Classify centerline
        no_change = False
        full_change = False
        div_change = False
        if (not in_p1) and (not in_p1):
            no_change = True
        elif in_p1 and in_p2:
            full_change = True
        else:
            div_change = True
            div_count += 1

        # Old line information
        n = line.GetNumberOfPoints()
        centerline_cell_array.InsertNextCell(n)
        radius_array_data = line.GetPointData().GetArray(radiusArrayName).GetTuple1

        # Iterate through centerline points
        k = 0
        for p in range(n):
            p_old = np.asarray(line.GetPoint(p))
            if no_change:
                dx = 0

            elif full_change:
                if id1 <= p <= id2:
                    p_new = np.asarray(new_cl.GetPoint(p - id1))
                    dx = p_new - p_old
                else:
                    dx = 0

            elif div_change:
                if p <= id1:
                    dx = 0
                elif p <= diverging_id[div_count][0]:
                    p_new = np.asarray(new_cl.GetPoint(p - id1))
                    dx = p_new - p_old
                else:
                    p_old = np.asarray(line.GetPoint(diverging_id[div_count][0]))
                    p_new = np.asarray(new_cl.GetPoint(diverging_id[div_count][0] - id1))
                    dx = p_new - p_old

            if not smooth_line:
                dx = - dx

            p_old = np.asarray(line.GetPoint(p))
            centerline_points.InsertNextPoint(p_old + dx)
            radius_array.SetTuple1(count, radius_array_data(p))
            centerline_cell_array.InsertCellPoint(count)
            count += 1

    centerline.SetPoints(centerline_points)
    centerline.SetLines(centerline_cell_array)
    centerline.GetPointData().AddArray(radius_array)

    return centerline


def read_command_line_curvature(input_path=None, output_path=None):
    """
    Read arguments from commandline and return all values in a dictionary.
    If input_path and output_path are not None, then do not parse command line, but
    only return default values.

    Args:
        input_path (str): Input file path, positional argument with default None.
        output_path (str): Output file path, positional argument with default None.
    """
    # Description of the script
    description = "Manipulates a selected part of a tubular geometry" + \
                  "by creating a sharper or smoother version of the input geometry, " + \
                  "determined by a smoothed version of the centerline representing the selected part."

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    # Add common arguments
    required = not (input_path is not None and output_path is not None)
    add_common_arguments(parser, required=required)

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

    # Set problem specific arguments
    parser.add_argument("-fl", "--smooth-factor-line", type=float, default=1.0,
                        help="Smoothing factor of centerline curve.")
    parser.add_argument("-it", "--iterations", type=int, default=100,
                        help="Smoothing iterations of centerline curve.")
    parser.add_argument("-sl", "--smooth-line", type=str2bool, default=True,
                        help="Smooths centerline if True, anti-smooths if False")

    # Parse paths to get default values
    if required:
        args = parser.parse_args()
    else:
        args = parser.parse_args(["-i" + input_path, "-o" + output_path])

    return dict(input_filepath=args.ifile, smooth=args.smooth,
                smooth_factor=args.smooth_factor, smooth_factor_line=args.smooth_factor_line,
                iterations=args.iterations, smooth_line=args.smooth_line, output_filepath=args.ofile,
                poly_ball_size=args.poly_ball_size, region_of_interest=args.region_of_interest,
                region_points=args.region_points, resampling_step=args.resampling_step,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point)


def main_curvature():
    manipulate_curvature(**read_command_line_curvature())


if __name__ == "__main__":
    manipulate_curvature(**read_command_line_curvature())

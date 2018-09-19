from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from common import *


def curvature_variations(input_filepath, smooth, smooth_factor_voro, smooth_factor_line, iterations,
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
        smooth_factor_voro (float): Smoothing factor used for Voronoi diagram smoothing.
        smooth_factor_line (float): Smoothing factor used for centerline smoothing.
        iterations (float): Smoothing iterations of centerline.
        region_of_interest (str): Method for setting the region of interest ['manuall' | 'commandline' | 'first_line']
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

    # Compute centerlines
    inlet, outlets = get_centers(surface, base_path)

    print("Compute centerlines and Voronoi diagram")
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=resampling_step,
                                                         smooth=False, base_path=base_path)
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines, base_path,
                                          smooth, smooth_factor_voro, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)
    # Get region of interest
    _, region_points = get_line_to_change(capped_surface, centerlines,
                                          region_of_interest, "variation", region_points, 0)
    region_points = [[region_points[3 * i], region_points[3 * i + 1], region_points[3 * i + 2]]
                     for i in range(len(region_points) // 3)]

    # Set and get clipping points, centerlines and diverging centerlines
    centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids = \
        find_region_of_interest_and_diverging_centerlines(centerlines, region_points)
    diverging_centerline_ispresent = diverging_centerlines is not None

    # Set region and extract relevant centerline
    locator = get_locator(extract_single_line(centerlines, 0))
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    p1 = centerlines.GetPoint(id1)
    p2 = centerlines.GetPoint(id2)
    centerline = extract_single_line(centerlines, 0, startID=id1, endID=id2)

    if diverging_centerline_ispresent:
        print("Clipping diverging centerlines")
        div_patch_cl = []
        for i, divline in enumerate(diverging_centerlines):
            clip_id = int((divline.GetNumberOfPoints() - diverging_ids[i]) * 0.2 + diverging_ids[i])
            patch_cl = clip_diverging_line(divline, centerlines.GetPoint(0), clip_id)
            div_patch_cl.append(extract_single_line(patch_cl, 1))

    # Clipp Voronoi diagram
    print("Clipping voronoi diagram")
    masked_voronoi = MaskVoronoiDiagram(voronoi, centerline)
    clipped_voronoi = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)

    # Extract voronoi diagram of diverging lines
    clipped_voronoi_div = []
    if diverging_centerline_ispresent:
        for div_line_ends in div_patch_cl:
            masked_voronoi_div = MaskVoronoiDiagram(voronoi, div_line_ends)
            clipped_voronoi_div.append(ExtractMaskedVoronoiPoints(voronoi, masked_voronoi_div))

    # Smooth relevant centerline
    smooth_carotid_siphon = vmtk_centerline_geometry(centerline, True, True,
                                                     factor=smooth_factor_line, iterations=iterations)
    write_polydata(smooth_carotid_siphon, centerline_smooth_path)

    # Get rest of artery
    p_start = centerline.GetPoint(0)
    p_end = centerline.GetPoint(centerline.GetNumberOfPoints() - 1)
    region = [p_start, p_end]
    region = np.asarray(region)
    region_vtk = vtk.vtkPoints()
    for point in region:
        region_vtk.InsertNextPoint(point)

    patch_cl = CreateParentArteryPatches(centerlines, region_vtk, siphon=True)
    end_lines = []
    for i in range(patch_cl.GetNumberOfCells()):
        tmp_line = extract_single_line(patch_cl, i)
        if tmp_line.GetNumberOfPoints() > 40:
            end_lines.append(tmp_line)

    centerlines_end = merge_data(end_lines)
    masked_voronoi_end = MaskVoronoiDiagram(voronoi, centerlines_end)
    clipped_voronoi_end = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi_end)

    # Update clipped Voronoi diagram
    print("Smooth / sharpen voronoi diagram")
    moved_clipped_voronoi = make_voronoi_smooth(clipped_voronoi, centerline, smooth_carotid_siphon, smooth_line)

    # Move diverging centerlines
    moved_clipped_voronoi_div = []
    for i, clipped_voronoi_div_i in enumerate(clipped_voronoi_div):
        moved_clipped_voronoi_div.append(make_voronoi_smooth(clipped_voronoi_div_i,
                                                             centerline,
                                                             smooth_carotid_siphon,
                                                             smooth_line, div=True,
                                                             div_point=region[i]))

    moved_clipped_voronoi_div.append(moved_clipped_voronoi)
    moved_clipped_voronoi_div.append(clipped_voronoi_end)

    # Combine Voronoi diagram
    new_voronoi = merge_data(moved_clipped_voronoi_div)

    # Create new surface
    print("Create new surface")
    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)
    new_centerlines = old_centerlines = merge_data([centerlines, diverging_centerlines])
    new_surface = prepare_surface_output(new_surface, surface,
                                         new_centerlines, output_filepath,
                                         test_merge=False, changed=True,
                                         old_centerline=old_centerlines)

    write_polydata(new_centerlines, new_centerlines_path)
    write_polydata(new_surface, output_filepath)


def make_voronoi_smooth(voronoi, old_cl, new_cl, smoothmode, div=False, div_point=None):
    """
    Move the voronoi diagram based on a smoothed
    version of the centerline.

    Args:
        voronoi (vtkPolyData): Voronoi diagram data set.
        old_cl (vtkPolyData): Unsmoothed centerline points.
        new_cl (vtkPolyData): Smoothed centerline points.
        smoothmode (bool): Determines if model becomes smoother or sharper.
        div (bool): True if centerline is a diverging line.
        div_point (ndarray): Diverging point along siphon.

    Returns:
        newDataSet (vtkPolyData): Manipulated voronoi diagram.
    """
    locator = get_locator(old_cl)
    N = voronoi.GetNumberOfPoints()
    newDataSet = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    # Define segments for transitioning
    IDend = old_cl.GetNumberOfPoints() - 1
    IDmidend = int(IDend * 0.9)
    IDstartmid = int(IDend * 0.1)
    IDmid = int(IDend * 0.2)
    # Iterate through voronoi points
    for i in range(N):
        if div:
            cl_id = locator.FindClosestPoint(div_point)
        else:
            cl_id = locator.FindClosestPoint(voronoi.GetPoint(i))
        p0 = np.asarray(old_cl.GetPoint(cl_id))
        p1 = np.asarray(new_cl.GetPoint(cl_id))
        if smoothmode:
            dx = p1 - p0
        else:
            dx = -(p1 - p0)

        # Smooth transition at inlet and at end of siphon
        if cl_id < IDstartmid:
            dx = 0
        elif IDstartmid <= cl_id < IDmid:
            dx = dx * (cl_id - IDstartmid) / float(IDmid - IDstartmid)
        elif cl_id > IDmidend:
            dx = dx * (IDend - cl_id) / float(IDend - IDmidend)

        dist = dx
        points.InsertNextPoint(np.asarray(voronoi.GetPoint(i)) + dist)
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)

    newDataSet.SetPoints(points)
    newDataSet.SetVerts(verts)
    newDataSet.GetPointData().AddArray(voronoi.GetPointData().GetArray(radiusArrayName))

    return newDataSet


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

    # Required arguments
    required.add_argument('-i', '--ifile', type=str, default=None,
                          help="Path to the surface model", required=True)
    required.add_argument("-o", "--ofile", type=str, default=None, required=True,
                          help="Relative path to the output surface. The default folder is" +
                               " the same as the input file, and a name with a combination of the" +
                               " parameters.")
    # Set region of interest:
    parser.add_argument("-r", "--region-of-interest", type=str, default="manuall",
                        choices=["manuall", "commandline", "first_line"],
                        help="The method for defining the region to be changed. There are" +
                             " three options: 'manuall', 'commandline', 'first_line'. In" +
                             " 'manuall' the user will be provided with a visualization of the" +
                             " input surface, and asked to provide an end and start point of the" +
                             " region of interest. Note that not all algorithms are robust over" +
                             " bifurcations. If 'commandline' is provided, then '--region-points'" +
                             " is expected to be provided. Finally, if 'first_line' is" +
                             " given, the method will chose the first line of the geometry.")
    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" +
                             " argument have to be given. The method expects two points" +
                             " which defines the start and end of the region of interest. ")
    # Optional arguments

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
    parser.add_argument('-s', '--smooth', type=bool, default=True,
                        help="Smooth the voronoi diagram, default is True")
    parser.add_argument('-f', '--smooth_factor', type=float, default=0.25,
                        help="If smooth option is true then each voronoi point" +
                             " that has a radius less then MISR*(1-smooth_factor) at" +
                             " the closest centerline point is removed.")
    parser.add_argument("-fl", "--smooth_factor_line", type=float, default=1.0,
                        help="Smoothing factor of centerline curve.")
    parser.add_argument("-it", "--iterations", type=int, default=100,
                        help="Smoothing iterations of centerline curve.")
    parser.add_argument("-sl", "--smooth_line", type=bool, default=True,
                        help="Smoothes centerline if True, anti-smoothes if False")
    parser.add_argument("-rs", "--resampling-step", type=float, default=0.1,
                        help="Resampling step for centerline resampling.")
    args = parser.parse_args()

    return dict(input_filepath=args.ifile, smooth=args.smooth,
                smooth_factor_voro=args.smooth_factor, smooth_factor_line=args.smooth_factor_line,
                iterations=args.iterations, smooth_line=args.smooth_line, output_filepath=args.ofile,
                poly_ball_size=args.poly_ball_size, region_of_interest=args.region_of_interest,
                region_points=args.region_points, resampling_step=args.resampling_step,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point)


if __name__ == "__main__":
    curvature_variations(**read_command_line())

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from scipy.ndimage.filters import gaussian_filter
from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from .common import *
#from argparse_common import *


def manipulate_area(input_filepath, method, smooth, smooth_factor, no_smooth,
                    no_smooth_point, region_of_interest, region_points, beta, ratio,
                    stenosis_length, percentage, output_filepath, poly_ball_size,
                    resampling_step):
    """
    Objective manipulation of area variation in
    patient-specific models of blood vessels.
    Manipulation is performed by looping over
    all Voronoi points relative to a selected centerline
    and change the corresponding radius.

    Args:
        input_filepath (str): Path to directory with cases.
        method (str): Type of manipulation of the centerline.
        smooth (bool): Determines Voronoi smoothing or not.
        smooth_factor (float): Smoothing factor used for voronoi diagram smoothing.
        no_smooth (bool): If True, define a region where the Voronoi diagram should not be smoothed.
        no_smooth_point (list): A flattened list to the 'end' points of the regions not to smooth.
        beta (float): Factor determining how much the geometry will differ.
        ratio (float): Target ratio, A_max / A_min. Beta is ignored (and estimated) if given.
        percentage (float): Percentage the area of the geometry / stenosis is increase/decreased.
        stenosis_length (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'first_line']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        poly_ball_size (list): Resolution of polyballs used to create surface.
        output_filepath (str): Path to output the manipulated surface.
        resampling_step (float): Resampling length for centerline resampling.
    """
    base_path = get_path_names(input_filepath)

    # Files paths
    voronoi_new_path = base_path + "_voronoi_manipulated.vtp"
    centerlines_path = base_path + "_centerline.vtp"
    centerline_area_spline_path = base_path + "_centerline_area_spline.vtp"
    centerline_area_spline_sections_path = base_path + "_centerline_area_sections.vtp"
    centerline_spline_path = base_path + "_centerline_spline.vtp"
    centerline_remaining_path = base_path + "_centerline_remaining.vtp"
    centerline_diverging_path = base_path + "_centerline_diverging.vtp"

    # Clean, triangulate, and capp/uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Create centerline and voronoi diagram
    inlet, outlets = get_centers(surface, base_path)
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=resampling_step,
                                                         smooth=False, base_path=base_path)

    # Smooth voronoi diagram
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)

    # Spline centerline and compute cross-sectional areas along line
    centerline_splined, centerline_remaining, \
        centerline_diverging, region_points = get_line_to_change(capped_surface, centerlines,
                                                                 region_of_interest, method,
                                                                 region_points, stenosis_length)
    write_polydata(centerline_splined, centerline_spline_path)
    write_polydata(centerline_remaining, centerline_remaining_path)
    if centerline_diverging is not None:
        write_polydata(merge_data(centerline_diverging), centerline_diverging_path)

    # Compute area
    centerline_area, centerline_area_sections = vmtk_compute_centerline_sections(surface,
                                                                                 centerline_splined)
    write_polydata(centerline_area, centerline_area_spline_path)
    write_polydata(centerline_area_sections, centerline_area_spline_sections_path)

    # Manipulate the voronoi diagram
    print("-- Change Voronoi diagram")
    centerline_regions = [centerline_splined, centerline_remaining]
    if centerline_diverging is not None:
        centerline_regions += centerline_diverging
    else:
        centerline_regions += [None]
    voronoi_regions = split_voronoi_with_centerlines(voronoi, centerline_regions)

    new_voronoi = change_area(voronoi_regions[0], centerline_area, method, beta, ratio, percentage,
                              region_of_interest, region_points, centerline_diverging,
                              voronoi_regions[2:])

    new_voronoi = merge_data([new_voronoi, voronoi_regions[1]])
    write_polydata(new_voronoi, voronoi_new_path)

    # Make new surface
    print("-- Create surface")
    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    print("-- Smoothing, clean, and check surface.")
    new_surface = prepare_surface_output(new_surface, surface, centerlines,
                                         output_filepath, test_merge=True)
    write_polydata(new_surface, output_filepath)


def get_factor(line_to_change, method, beta, ratio, percentage, region_of_interest,
               region_points):
    """
    Compute the factor determining
    the change in radius, used to
    manipulate the Voronoi diagram.

    Args:
        line_to_change (vtkPolyData): Centerline representing area of interest.
        method (str): Type of manipulation of the centerline.
        beta (float): Factor deciding how area will change. Ignored if ratio is given.
        ratio (float): Desired ratio between min and max cross-sectional area.
        percentage (float): Desired increase/decrease in cross-sectional area.
        region_of_interest (str): Method for setting the region of interest.
        region_points (list): List of points for the stenosis.

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Array to change the radius
    area = get_array("CenterlineSectionArea", line_to_change)

    # Safety smoothing, section area does not always work perfectly
    for i in range(2):
        area = gaussian_filter(area, 5)
    mean_area = np.mean(area)

    # Linear transition first and last 10 % for some combinations of method an region_of_interest
    if region_of_interest in ["manual", "commandline"] and method in ["area", "variation"]:
        k = int(round(area.shape[0] * 0.10, 0))
        linear = area.shape[0] - k * 2
    else:
        k = 0
        linear = area.shape[0]

    # Transition
    trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(linear).tolist() +
                       np.linspace(0, 1, k).tolist())

    # Only smooth end with first_line
    if region_of_interest == "first_line":
        k = int(round(area.shape[0] * 0.10, 0))
        linear = area.shape[0] - k
        trans = np.asarray(np.zeros(linear).tolist() + np.linspace(0, 1, k).tolist())

    # Get factor
    if method == "variation":
        if ratio is not None:
            # Initial guess
            R_old = area.max() / area.min()
            beta = 0.5 * (math.log(ratio) / math.log(R_old) - 1)
            factor_ = (area / mean_area) ** beta
            factor = factor_[:, 0] * (1 - trans) + trans
        else:
            factor_ = ((area / mean_area) ** beta)[:, 0]
            factor = factor_ * (1 - trans) + trans

    elif method == "stenosis":
        if len(region_points) == 3:
            t = np.linspace(0, np.pi, line_to_change.GetNumberOfPoints())
            factor = (1 - np.sin(t) * percentage * 0.01).tolist()
            factor = factor * (1 - trans) + trans

        elif len(region_points) == 6:
            length = get_curvilinear_coordinate(line_to_change)
            factor = (area[0] + (area[-1] - area[0]) * (length / length.max())) / area[:, 0]
            factor = np.sqrt(factor)
            factor = factor * (1 - trans) + trans

    # Increase or decrease overall area by a percentage
    elif method == "area":
        factor_ = np.ones(len(trans)) * (1 + percentage * 0.01)
        factor = factor_ * (1 - trans) + trans

    return factor


def change_area(voronoi, line_to_change, method, beta, ratio, percentage,
                region_of_interest, region_points, diverging_centerline,
                diverging_voronoi):
    """
    Change the cross-sectional area of an input
    voronoi diagram along the corresponding area
    represented by a centerline.

    Args:
        voronoi (vtkPolyData): Voronoi diagram.
        line_to_change (vtkPolyData): Centerline representing area of interest.
        method (str): Type of manipulation of the centerline.
        beta (float): Factor deciding how area will change. Ignored if ratio is given.
        ratio (float): Desired ratio between min and max cross-sectional area.
        percentage (float): Percentage the area of the geometry / stenosis is increase/decreased.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'first_line']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        diverging_centerline (vtkPolyData): Polydata containing diverging centerlines along region of interest.
        diverging_voronoi (vtkPolyData): Voronoi diagram diverging off region of interest.

    Returns:
        new_voronoi (vtkPolyData): Manipulated Voronoi diagram.
    """
    # Get factor
    factor = get_factor(line_to_change, method, beta, ratio, percentage,
                        region_of_interest, region_points)

    # Locator to find closest point on centerline
    locator = get_locator(line_to_change)

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

        tmp_id = locator.FindClosestPoint(point)

        v1 = np.asarray(line_to_change.GetPoint(tmp_id)) - np.asarray(point)
        v2 = v1 * (1 - factor[tmp_id])
        point = (np.asarray(point) + v2).tolist()

        # Change radius
        point_radius = point_radius * factor[tmp_id]

        voronoi_points.InsertNextPoint(point)
        radius_array.SetTuple1(i, point_radius)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i)

    # Offset Voronoi diagram along "diverging" centerlines
    if diverging_centerline is not None:
        count = i
        for j in range(len(diverging_voronoi)):
            # Get offset for Voronoi diagram
            for i in range(len(diverging_centerline.GetNumberOfPoints())):
                point = diverging_centerline.GetPoint(i)
                tmp_id = locator.FindClosestPoint(point)
                dist = distance(diverging_centerline.GetPoint(i), line_to_change.GetPoint(tmp_id))
                if dist > line_to_change.GetPointData().GetArray(radiusArrayName).GetTuple1(tmp_id):
                    break

            v1 = np.array(line_to_change.GetPoint(tmp_id)) - np.array(point)
            v2 = v1 * (1 - factor[tmp_id])

            # Offset Voronoi diagram
            for i in range(diverging_voronoi):
                point = diverging_voronoi[j].GetPoint(i)
                point_radius = diverging_voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
                point = (np.asarray(point) + v2).tolist()

                voronoi_points.InsertNextPoint(point)
                radius_array.SetTuple1(count, point_radius)
                cell_array.InsertNextCell(1)
                cell_array.InsertCellPoint(count)
                count += 1

    new_voronoi.SetPoints(voronoi_points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)

    return new_voronoi


def read_command_line():
    """
    Read arguments from commandline
    """
    # Description of the script
    description = "Manipulates the area of a tubular geometry. The script changes the area" + \
                  " in three different ways:" + \
                  "\n1) Increase or decrease the area variation along the region of" + \
                  " interest. (variation)\n2) Create or remove a local narrowing. (stenosis)" + \
                  "\n3) Inflate or deflate the entire region of interest. (area)"
    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    # Add common arguments
    add_common_arguments(parser)

    # Mode / method for manipulation
    parser.add_argument("-m", "--method", type=str, default="variation",
                        choices=["variation", "stenosis", "area"],
                        help="Methods for manipulating the area in the region of interest:" +
                             "\n1) 'variation' will increase or decrease the changes in area" +
                             " along the centerline of the region of interest." +
                             "\n2) 'stenosis' will create or remove a local narrowing of the" +
                             " surface. If two points is provided, the area between these" +
                             " two points will be linearly interpolated to remove the narrowing." +
                             " If only one point is provided it is assumed to be the center of" +
                             " the stenosis. The new stenosis will have a sin shape, however, any" +
                             " other shape may be easily implemented." +
                             "\n3) 'area' will inflate or deflate the area in the region of" +
                             " interest.")

    # Set region of interest
    parser.add_argument("-r", "--region-of-interest", type=str, default="manual",
                        choices=["manual", "commandline", "first_line"],
                        help="The method for defining the region to be changed. There are" +
                             " three options: 'manual', 'commandline', 'first_line'. In" +
                             " 'manual' the user will be provided with a visualization of the" +
                             " input surface, and asked to provide an end and start point of the" +
                             " region of interest. Note that not all algorithms are robust over" +
                             " bifurcations. If 'commandline' is provided, then '--region-points'" +
                             " is expected to be provided. Finally, if 'first_line' is given, the" +
                             " line from the inlet (largest opening) to the first bifurcation" +
                             " will be altered, not that method='stenosis' can not be used" +
                             " with 'first_line'.")
    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" +
                             " argument have to be given. The method expects two points" +
                             " which defines the start and end of the region of interest. If" +
                             " 'method' is set to stenosis, then one point can be provided as well," +
                             " which is assumed to be the center of a new stenosis." +
                             " Example providing the points (1, 5, -1) and (2, -4, 3):" +
                             " --stenosis-points 1 5 -1 2 -4 3")

    # "Variation" arguments
    parser.add_argument('--beta', type=float, default=0.5,
                        help="For method=variation: The new voronoi diagram is computed as" +
                             " (A/A_mean)**beta*r_old, over the respective area. If beta <" +
                             " 0 the geometry will have less area variation, and" +
                             " if beta > 0, the variations in area will increase")
    parser.add_argument("--ratio", type=float, default=None,
                        help="For method=variation: " +
                             " Target ratio of A_max/A_min, when this is given beta will be ignored" +
                             " and instead approximated to obtain the target ratio")

    # "Stenosis" argument
    parser.add_argument("--size", type=float, default=2.0, metavar="length",
                        help="For method=stenosis: The length of the area " +
                             " affected by a stenosis relative to the minimal inscribed" +
                             " sphere radius of the selected point. Default is 2.0.")

    # "area" / "stenosis" argument
    parser.add_argument("--percentage", type=float, default=50.0,
                        help="Percentage the" +
                             " area of the geometry is increase/decreased overall or only" +
                             " stenosis")

    # Parse
    args = parser.parse_args()

    if args.method == "stenosis" and args.region_of_interest == "first_line":
        raise ValueError("Can not set region of interest to 'first_line' when creating or" +
                         " removing a stenosis")

    if args.method == "variation" and args.ratio is not None and args.beta != 0.5:
        print("WARNING: The beta value you provided will be ignored, using ratio instead.")

    if args.region_points is not None:
        if len(args.region_points) % 3 != 0 or len(args.region_points) > 6:
            raise ValueError("ERROR: Please provide region point(s) as a multiple of 3, and maximum" +
                             " two points.")

    if args.no_smooth_point is not None and len(args.no_smooth_point):
        if len(args.no_smooth_point) % 3 != 0:
            raise ValueError("ERROR: Please provide the no smooth point(s) as a multiple" +
                             " of 3.")

    return dict(input_filepath=args.ifile, method=args.method, smooth=args.smooth,
                smooth_factor=args.smooth_factor, beta=args.beta,
                region_of_interest=args.region_of_interest,
                region_points=args.region_points, ratio=args.ratio,
                stenosis_length=args.size,
                percentage=args.percentage, output_filepath=args.ofile,
                poly_ball_size=args.poly_ball_size, no_smooth=args.no_smooth,
                no_smooth_point=args.no_smooth_point, resampling_step=args.resampling_step)


def main_area():
    manipulate_area(**read_command_line())

if __name__ == '__main__':
    manipulate_area(**read_command_line())

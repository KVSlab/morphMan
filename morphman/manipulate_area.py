##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from morphman.common.argparse_common import *
from morphman.common.surface_operations import *


def manipulate_area(input_filepath, method, smooth, smooth_factor, no_smooth,
                    no_smooth_point, region_of_interest, region_points, beta, ratio,
                    size, percentage, output_filepath, poly_ball_size,
                    resampling_step, angle_asymmetric):
    """
    Objective manipulation of area variation in
    patient-specific models of blood vessels.
    Manipulation is performed by looping over
    all Voronoi points relative to a selected centerline
    and change the corresponding radius.

    Args:
        angle_asymmetric (float): Angle defining the orientation of the asymmetry for a stenosis / bulge
        input_filepath (str): Path to directory with cases.
        method (str): Type of manipulation of the centerline.
        smooth (bool): Determines Voronoi smoothing or not.
        smooth_factor (float): Smoothing factor used for voronoi diagram smoothing.
        no_smooth (bool): If True, define a region where the Voronoi diagram should not be smoothed.
        no_smooth_point (list): A flattened list to the 'end' points of the regions not to smooth.
        beta (float): Factor determining how much the geometry will differ.
        ratio (float): Target ratio, A_max / A_min. Beta is ignored (and estimated) if given.
        percentage (float): Percentage the area of the geometry / stenosis is increase/decreased.
        size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'first_line']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        poly_ball_size (list): Resolution of polyballs used to create surface.
        output_filepath (str): Path to output the manipulated surface.
        resampling_step (float): Resampling length for centerline resampling.
    """
    base_path = get_path_names(input_filepath)

    # Files paths
    centerlines_path = base_path + "_centerline.vtp"
    centerline_area_spline_path = base_path + "_centerline_area_spline.vtp"
    centerline_area_spline_sections_path = base_path + "_centerline_area_sections.vtp"
    centerline_spline_path = base_path + "_centerline_spline.vtp"
    centerline_remaining_path = base_path + "_centerline_remaining.vtp"
    centerline_diverging_path = base_path + "_centerline_diverging.vtp"

    voronoi_new_path = base_path + "_voronoi_manipulated.vtp"
    voronoi_roi_path = base_path + "_voronoi_region_of_intrest.vtp"
    voronoi_rest_path = base_path + "_voronoi_rest.vtp"
    voronoi_div_path = base_path + "_voronoi_div{:d}.vtp"

    # Clean, triangulate, and capp/uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Create centerline and voronoi diagram
    inlet, outlets = get_inlet_and_outlet_centers(surface, base_path)
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=resampling_step,
                                                         smooth=False, base_path=base_path)

    # Smooth voronoi diagram
    if smooth:
        voronoi = prepare_voronoi_diagram(capped_surface, centerlines, base_path, smooth,
                                          smooth_factor, no_smooth, no_smooth_point,
                                          voronoi, pole_ids, resampling_step)

    # Spline centerline and compute cross-sectional areas along line
    centerline_splined, centerline_remaining, centerline_diverging, region_points, diverging_ids = get_line_to_change(
        capped_surface, centerlines,
        region_of_interest, method,
        region_points, size)
    write_polydata(centerline_splined, centerline_spline_path)
    write_polydata(centerline_remaining, centerline_remaining_path)
    if centerline_diverging is not None:
        write_polydata(vtk_merge_polydata(centerline_diverging), centerline_diverging_path)

    # Split the Voronoi diagram
    print("-- Changing Voronoi diagram")
    centerline_regions = [centerline_splined, centerline_remaining]
    if centerline_diverging is not None:
        for i, div_cl in enumerate(centerline_diverging):
            centerline_regions += [extract_single_line(div_cl, 0, start_id=diverging_ids[i])]
    voronoi_regions = get_split_voronoi_diagram(voronoi, centerline_regions)

    # Write the seperate segments
    write_polydata(voronoi_regions[0], voronoi_roi_path)
    write_polydata(voronoi_regions[1], voronoi_rest_path)
    for i in range(2, len(voronoi_regions)):
        write_polydata(voronoi_regions[i], voronoi_div_path.format(i - 1))

    # Compute area, and acount for diverging branches.
    if centerline_diverging is not None or smooth:
        surface_area = create_new_surface(voronoi_regions[0],
                                          poly_ball_size=poly_ball_size)
    else:
        surface_area = surface
    centerline_area, centerline_area_sections = vmtk_compute_centerline_sections(surface_area,
                                                                                 centerline_splined)
    write_polydata(centerline_area, centerline_area_spline_path)
    write_polydata(centerline_area_sections, centerline_area_spline_sections_path)

    # Manipulate the Voronoi diagram
    factor = get_factor(centerline_area, method, beta, ratio, percentage,
                        region_of_interest)
    new_voronoi, new_centerlines = change_area(voronoi_regions[0], factor, centerline_area,
                                               centerline_diverging, voronoi_regions[2:],
                                               surface_area, centerlines, angle_asymmetric)
    new_voronoi = vtk_merge_polydata([new_voronoi, voronoi_regions[1]])
    write_polydata(new_voronoi, voronoi_new_path)

    # Make new surface
    print("-- Creating surface")
    new_surface = create_new_surface(new_voronoi, poly_ball_size=poly_ball_size)

    print("-- Smoothing, cleaning, and checking surface.")
    new_surface = prepare_output_surface(new_surface, surface, centerlines,
                                         output_filepath, test_merge=True)

    print("\n-- Writing new surface to {}.".format(output_filepath))
    write_polydata(new_surface, output_filepath)


def get_factor(line_to_change, method, beta, ratio, percentage, region_of_interest):
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

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Array to change the radius
    area = get_point_data_array("CenterlineSectionArea", line_to_change)
    mean_area = np.mean(area)

    # Linear transition first and last 10 % for some combinations of method an region_of_interest
    if region_of_interest in ["manual", "commandline"] and method in ["area", "variation"]:
        k = int(round(area.shape[0] * 0.10, 0))
        linear = area.shape[0] - k * 2
    else:
        k = 0
        linear = area.shape[0]

    # Transition
    trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(linear).tolist() + np.linspace(0, 1, k).tolist())

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

    # Create or remove stenosis
    elif method == "stenosis":
        t = np.linspace(0, np.pi, line_to_change.GetNumberOfPoints())
        factor = (1 - np.sin(t) * percentage * 0.01).tolist()
        factor = factor * (1 - trans) + trans

    elif method == "linear":
        # Note: No need for a transition since the area should not change in the start and end
        length = get_curvilinear_coordinate(line_to_change)
        factor = (area[0] + (area[-1] - area[0]) * (length / length.max())) / area[:, 0]
        factor = np.sqrt(factor)

    # Create a fusiform aneurysm or bulge
    elif method == "bulge":
        t = np.linspace(0, np.pi, line_to_change.GetNumberOfPoints())
        factor = (1 + np.sin(t) * percentage * 0.01).tolist()
        factor = factor * (1 - trans) + trans

    # Increase or decrease overall area by a percentage
    elif method == "area":
        factor_ = np.ones(len(trans)) * (1 + percentage * 0.01)
        factor = factor_ * (1 - trans) + trans

    return factor


def change_area(voronoi, factor, line_to_change, diverging_centerline, diverging_voronoi, surface_area, centerlines,
                angle_asymmetric):
    """
    Change the cross-sectional area of an input
    voronoi diagram along the corresponding area
    represented by a centerline.

    Args:
        angle_asymmetric (float): Angle defining the orientation of the asymmetry for a stenosis / bulge
        voronoi (vtkPolyData): Voronoi diagram.
        factor (ndarray): An array with a factor for changing each point along the centerline
        line_to_change (vtkPolyData): Centerline representing area of interest.
        diverging_centerline (list): List of polydata containing diverging centerlines along region of interest.
        diverging_voronoi (list): List of Voronoi diagram diverging off region of interest.
        surface_area (vtkPolyData): The surface used to compute the area.
        centerlines (vtkPolyData): Centerlines of the full model.

    Returns:
        new_voronoi (vtkPolyData): Manipulated Voronoi diagram.
    """
    # Locator to find closest point on centerline
    locator = get_vtk_point_locator(line_to_change)

    # Voronoi diagram
    n = voronoi.GetNumberOfPoints()
    if len(diverging_voronoi) > 0:
        m = n + sum([div_voro.GetNumberOfPoints() for div_voro in diverging_voronoi])
    else:
        m = n

    new_voronoi = vtk.vtkPolyData()
    voronoi_points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array = get_vtk_array(radiusArrayName, 1, m)
    N = line_to_change.GetNumberOfPoints()

    # Iterate through Voronoi diagram and manipulate
    frenet_normals_array = get_point_data_array("FrenetNormal", line_to_change, k=3)
    frenet_tangents_array = get_point_data_array("FrenetTangent", line_to_change, k=3)
    point_radius_array = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1
    tol = get_centerline_tolerance(centerlines)
    for i in range(n):
        id_list = vtk.vtkIdList()
        point = voronoi.GetPoint(i)

        # Find two closest points on centerline
        locator.FindClosestNPoints(2, point, id_list)
        tmp_id1, tmp_id2 = id_list.GetId(0), id_list.GetId(1)

        # Define points and vectors
        A = np.asarray(line_to_change.GetPoint(tmp_id1))
        P = np.asarray(point)
        B = np.asarray(line_to_change.GetPoint(tmp_id2))
        AB = B - A
        AP = P - A

        # If point is "beyond" start or end of centerline
        sign = 1
        if tmp_id1 in [0, N]:
            tetha = get_angle(AP, AB)
            if tetha > np.pi / 2.0:
                sign = -1

        # Get direction (delta_p) to move the point
        AB_length = la.norm(AB)
        if AB_length > tol * 10:
            h = np.linalg.norm(np.cross(AP, AB)) / AB_length  # shortest distance between point and line
            D = A + (AB / AB_length) * np.sqrt(np.linalg.norm(AP) ** 2 - h ** 2)
        else:
            D = A

        delta_p = D - P

        # Update factor value based on midpoint
        if sign > 0:
            factor_ = update_factor(A, AB_length, B, D, factor, tmp_id1, tmp_id2)
        else:
            factor_ = factor[tmp_id1]

        if angle_asymmetric is not None:
            v = get_asymmetric_displacement(A, angle_asymmetric, factor_, frenet_normals_array,
                                            frenet_tangents_array, locator, point)
        else:
            v = delta_p * (1 - factor_)

        # Change radius
        point_radius = point_radius_array(i) * factor_
        voronoi_points.InsertNextPoint((P + v).tolist())
        radius_array.SetTuple1(i, point_radius)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i)

    # Offset Voronoi diagram along "diverging" centerlines
    if diverging_centerline is not None:
        count = i + 1
        loc_surf = get_vtk_point_locator(surface_area)
        loc_cl = get_vtk_point_locator(centerlines)

        # 'Copy' old centerlines
        new_centerlines = vtk.vtkPolyData()
        new_centerlines.SetPoints(centerlines.GetPoints())
        new_centerlines.SetVerts(centerlines.GetVerts())
        new_centerlines.SetLines(centerlines.GetLines())
        new_centerlines.GetPointData().AddArray(centerlines.GetPointData().GetArray(radiusArrayName))
        points = new_centerlines.GetPoints()
        for j in range(len(diverging_voronoi)):
            # Get closest point on centerline to surface
            min_dist = 1e10
            for i in range(diverging_centerline[j].GetNumberOfPoints()):
                point = diverging_centerline[j].GetPoint(i)
                dist = get_distance(surface_area.GetPoint(loc_surf.FindClosestPoint(point)), point)
                if dist < min_dist:
                    min_point = point
                    min_dist = dist
                    div_id = i

            # Get offset for Voronoi diagram and centerline
            tmp_id1, tmp_id2 = id_list.GetId(0), id_list.GetId(1)

            A = np.asarray(line_to_change.GetPoint(tmp_id1))
            B = np.asarray(min_point)
            C = np.asarray(line_to_change.GetPoint(tmp_id2))
            BA = A - B
            BC = C - B
            AC = C - A
            AC_length = np.linalg.norm(AC)

            if AC_length != 0:
                h = np.linalg.norm(np.cross(BA, BC)) / AC_length  # shortest distance between point and line
                D = A + (AC / AC_length) * np.sqrt(np.linalg.norm(BA) ** 2 - h ** 2)
            else:
                D = A

            delta_p = D - B
            factor_ = factor[tmp_id1] * (1 - np.linalg.norm(D) / AC_length) \
                      + factor[tmp_id2] * (np.linalg.norm(D) / AC_length)
            v = delta_p * (1 - factor_)

            # Offset centerline
            # This method might not work if diverging centerline has a daughter branch
            for k in range(div_id, diverging_centerline[j].GetNumberOfPoints()):
                point = diverging_centerline[j].GetPoint(k)
                cl_id = loc_cl.FindClosestPoint(point)
                points.SetPoint(cl_id, (np.array(point) + v).tolist())

            # Offset Voronoi diagram
            point_radius_array = diverging_voronoi[j].GetPointData().GetArray(radiusArrayName).GetTuple1
            for i in range(diverging_voronoi[j].GetNumberOfPoints()):
                point = diverging_voronoi[j].GetPoint(i)
                point = (np.asarray(point) + v).tolist()

                voronoi_points.InsertNextPoint(point)
                radius_array.SetTuple1(count, point_radius_array(i))
                cell_array.InsertNextCell(1)
                cell_array.InsertCellPoint(count)
                count += 1

        new_centerlines.SetPoints(points)

    else:
        new_centerlines = centerlines

    new_voronoi.SetPoints(voronoi_points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)

    return new_voronoi, new_centerlines


def get_asymmetric_displacement(A, angle_asymmetric, factor_, frenet_normals_array, frenet_tangents_array, locator,
                                point):
    """
    Compute translation adding asymmetric effect to stenosis / bulge

    Args:
        A (ndarray): Initial position on centerline
        angle_asymmetric (Angle to rotate around Frenet tangent vector:
        factor_ (float): Scaling factor for area manipulation
        frenet_normals_array (list): List of Frenet normal vectors
        frenet_tangents_array (list): List of Frenet tangent vectors
        locator (vtkPointLocator): Centerline locator
        point (vtkPoint): Current Voronoi point

    Returns:
        ndarray: Translation vector for asymmetric effect
    """
    cl_id = locator.FindClosestPoint(point)
    frenet_normal = frenet_normals_array[cl_id]
    frenet_normal_vector = A + frenet_normal
    frenet_normal_vector /= la.norm(frenet_normal_vector)
    frenet_tangent = frenet_tangents_array[cl_id]
    R = get_rotation_matrix(frenet_tangent, angle_asymmetric)
    v = np.dot(R, frenet_normal_vector) * (1 - factor_) * 1.5

    return v


def update_factor(A, AB_length, B, P_mid, factor, tmp_id1, tmp_id2):
    """
    Update values in Factor based on midpoint between
    A and B.

    Args:
        P_mid (ndarray): Midpoint
        factor (list): List of scaling factors for area variation
        tmp_id1 (int): Id at first point
        tmp_id2 (int): Id at second point
        A (ndarray): First point
        B (ndarray): Second point
        AB_length (float): Length of line AB

    Returns:
        list: List of updated scaling factors
    """
    AP_mid_length = la.norm(P_mid - A)
    BP_mid_length = la.norm(P_mid - B)
    factor_ = (factor[tmp_id1] * AP_mid_length + factor[tmp_id2] * BP_mid_length) / AB_length

    return factor_


def read_command_line_area(input_path=None, output_path=None):
    """
    Read arguments from commandline and return all values in a dictionary.
    If input_path and output_path are not None, then do not parse command line, but
    only return default values.

    Args:
        input_path (str): Input file path, positional argument with default None.
        output_path (str): Output file path, positional argument with default None.
    """
    # Description of the script
    description = "Manipulates the area of a tubular geometry. The script changes the area" + \
                  " in three different ways:" + \
                  "\n1) Increase or decrease the area variation along the region of" + \
                  " interest. (variation)\n2) Create a local narrowing (stenosis) or widening (bulge)." + \
                  "\n3) Inflate or deflate the entire region of interest. (area)" + \
                  "\n4) A linear change in area between two points, for instance to remove a stenosis."
    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    # Add common arguments
    required = not (input_path is not None and output_path is not None)
    add_common_arguments(parser, required=required)

    # Mode / method for manipulation
    parser.add_argument("-m", "--method", type=str, default="variation",
                        choices=["variation", "stenosis", "area", "bulge", "linear"],
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
                        help="For method=[stenosis, bulge]: The length of the area " +
                             " affected by a stenosis relative to the minimal inscribed" +
                             " sphere radius of the selected point. Default is 2.0.")

    # "area" / "stenosis" argument
    parser.add_argument("--percentage", type=float, default=50.0,
                        help="Percentage the" +
                             " area of the geometry is increase/decreased overall or only" +
                             " in stenosis / bulge.")

    # Arguments for rotation of asymmetric area variation
    parser.add_argument('-as', '--angle-asymmetric', type=float, default=None,
                        help="Angle in degrees, defining asymmetric manipulation. " +
                             "Introduces asymmetric stenosis / bulges by applying an angle-dependent profile to " +
                             "area variation factor. Intended for 'stenosis' and 'bulge' methods, " +
                             "experimental for other methods.", metavar="angle-asymmetric")

    # Parse paths to get default values
    if required:
        args = parser.parse_args()
    else:
        args = parser.parse_args(["-i" + input_path, "-o" + output_path])

    if args.method in ["stenosis", "bulge", "linear"] and args.region_of_interest == "first_line":
        raise ValueError("Can not set region of interest to 'first_line' for 'stenosis'," + \
                         " 'bulge', or 'linear'")

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
    if args.angle_asymmetric is not None:
        angle_radians = args.angle_asymmetric * np.pi / 180.0  # Convert from deg to rad
    else:
        angle_radians = None

    return dict(input_filepath=args.ifile, method=args.method, smooth=args.smooth,
                smooth_factor=args.smooth_factor, beta=args.beta,
                region_of_interest=args.region_of_interest, angle_asymmetric=angle_radians,
                region_points=args.region_points, ratio=args.ratio,
                size=args.size, percentage=args.percentage, output_filepath=args.ofile,
                poly_ball_size=args.poly_ball_size, no_smooth=args.no_smooth,
                no_smooth_point=args.no_smooth_point, resampling_step=args.resampling_step)


def main_area():
    manipulate_area(**read_command_line_area())


if __name__ == '__main__':
    manipulate_area(**read_command_line_area())

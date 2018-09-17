from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelextrema
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from os import path, listdir

# Local import
from common import *

def area_variations(input_filepath, method, smooth, smooth_factor, no_smooth,
                    no_smooth_point, region_of_interest, region_points, beta, ratio,
                    stenosis_length, percentage, output_filepath, poly_ball_size):
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
        no_smooth_pooint (list): A flattend list to the 'end' points of the regions not to smooth.
        beta (float): Factor determining how much the geometry will differ.
        ratio (float): Target ratio, A_max / A_min. Beta is ignored (and estimated) if given.
        percentage (float): Percentage the area of the geometry / stenosis is increase/decreased.
        stenosis_length (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        region_of_interest (str): Method for setting the region of interest ['manuall' | 'commandline' | 'first_line']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
    """
    base_path = get_path_names(input_filepath)

    # Files paths
    surface_capped_path = base_path + "_capped.vtp"
    voronoi_new_path = base_path + "_voronoi_manipulated.vtp"
    centerlines_path = base_path + "_centerline.vtp"
    centerline_area_path = base_path + "_centerline_area.vtp"
    centerline_area_spline_path = base_path + "_centerline_area_spline.vtp"
    centerline_area_spline_sections_path = base_path + "_centerline_area_sections.vtp"
    centerline_spline_path = base_path + "_centerline_spline.vtp"

    # Clean, triangulate, and capp/uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Create centerline and voronoi diagram
    inlet, outlets = get_centers(surface, base_path)
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=0.1,
                                                         smooth=False, base_path=base_path)

    # Smooth voronoi diagram
    if smooth:
        voronoi = prepare_voronoi_diagram(surface, capped_surface, centerlines, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)

    # Spline centerline and compute cross-sectional areas along line
    centerline_splined, region_points = get_line_to_change(capped_surface, centerlines,
                                                           region_of_interest, method,
                                                           region_points, stenosis_length)
    write_polydata(centerline_splined, centerline_spline_path)
    centerline_area, centerline_area_sections = vmtk_compute_centerline_sections(surface,
                                                                                    centerline_splined)
    write_polydata(centerline_area, centerline_area_spline_path)
    write_polydata(centerline_area_sections, centerline_area_spline_sections_path)

    # Manipulate the voronoi diagram
    print("Change Voronoi diagram")
    newvoronoi = change_area(voronoi, centerline_area, method, beta, ratio, percentage,
                             region_of_interest, region_points, stenosis_length, surface)
    write_polydata(newvoronoi, voronoi_new_path)

    # Make new surface
    print("Create surface")
    new_surface = create_new_surface(newvoronoi, poly_ball_size=poly_ball_size)

    print("Cleaning surface for output")
    new_surface = prepare_surface_output(new_surface, surface, centerlines,
                                         output_filepath, test_merge=True)

    print("Write surface to: {}".format(output_filepath))
    write_polydata(new_surface, output_filepath)


def get_line_to_change(surface, centerline, region_of_interest, method, region_points,
                       stenosis_length):
    """
    Extract and spline part of centerline
    within the geometry where
    area variations will be performed.

    Args:
        centerline (vtkPolyData): Centerline in geometry.
        region_of_interest (str): Method for setting the region of interest ['manuall' | 'commandline' | 'first_line']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint.

    Returns:
        line_to_change (vtkPolyData): Part of centerline.
    """
    if region_of_interest == "first_line":
        tol = get_tolerance(centerline)
        line2 = extract_single_line(centerline, 0)
        numberOfPoints2 = line2.GetNumberOfPoints()

        # Iterate through lines and find diverging point
        n = 2
        pointIDs = []
        for j in range(1, n):
            line1 = extract_single_line(centerline, j)
            numberOfPoints1 = line1.GetNumberOfPoints()

            N = min(numberOfPoints1, numberOfPoints2)
            for i in range(N):
                point1 = line1.GetPoints().GetPoint(i)
                point2 = line2.GetPoints().GetPoint(i)
                if distance(point1, point2) > tol:
                    pointID = i
                    break
            pointIDs.append(pointID)

        startID = 0
        endID = min(pointIDs)
        cl_id = 0

        region_points = []

    elif region_of_interest == "commandline" or region_of_interest == "manuall":
        # Get points from the user
        if region_of_interest == "manuall":
            stenosis_point_id = vtk.vtkIdList()
            first = True
            while stenosis_point_id.GetNumberOfIds() not in [1, 2]:
                if not first:
                    print("Please provide only one or two points, try again")

                # Select point on surface
                seed_selector = vmtkPickPointSeedSelector()
                seed_selector.SetSurface(surface)
                if method == "variation" or method == "area":
                    seed_selector.text = "Press space to select the start and endpoint of the" + \
                                        " region of interest, 'u' to undo.\n"
                elif method == "stenosis":
                    seed_selector.text = "Press space to select, the center of a new" + \
                                         " stenosis (one point),\nOR place two points on each side" + \
                                         " of an existing stenosis to remove it, \'u\' to undo."
                seed_selector.Execute()
                stenosis_point_id = seed_selector.GetTargetSeedIds()
                first = True
            region_points = []
            for i in range(stenosis_point_id.GetNumberOfIds()):
                region_points += surface.GetPoint(stenosis_point_id.GetId(i))

        # Get locator
        locator = get_locator(centerline)

        if len(region_points) == 3:
            # Project point onto centerline
            point1 = region_points
            point1 = centerline.GetPoint(locator.FindClosestPoint(point1))

            # Get relevant line
            tol = get_tolerance(centerline)
            cl_id = -1
            dist = 1e10
            while dist > tol / 10:
                cl_id += 1
                line = extract_single_line(centerline, cl_id)
                tmp_loc = get_locator(line)
                tmp_id = tmp_loc.FindClosestPoint(point1)
                dist = distance(point1, line.GetPoint(tmp_id))

            # Get length of stenosis
            misr = get_array(radiusArrayName, line)
            length = stenosis_length * misr[tmp_loc.FindClosestPoint(point1)]

            # Get ids of start and stop
            centerline_length = get_curvilinear_coordinate(line)
            center = centerline_length[tmp_id]
            region_of_interest_id = (center - length <= centerline_length) \
                                 * (centerline_length <= center + length)
            startID = np.argmax(region_of_interest_id)
            endID = region_of_interest_id.shape[0] - 1 - np.argmax(region_of_interest_id[::-1])

        else:
            point1 = region_points[:3]
            point2 = region_points[3:]
            point1 = centerline.GetPoint(locator.FindClosestPoint(point1))
            point2 = centerline.GetPoint(locator.FindClosestPoint(point2))

            distance1 = []
            distance2 = []
            ids1 = []
            ids2 = []
            for i in range(centerline.GetNumberOfLines()):
                line = extract_single_line(centerline, i)
                tmp_loc = get_locator(line)
                ids1.append(tmp_loc.FindClosestPoint(point1))
                ids2.append(tmp_loc.FindClosestPoint(point2))
                cl_point1 = line.GetPoint(ids1[-1])
                cl_point2 = line.GetPoint(ids[-1])
                distance1.append(distance(point1, cl_point1))
                distance2.append(distance(point2, cl_point2))

            tol = get_tolerance(centerlines) / 10
            total_distance = (np.array(distance1) < tol) * (np.array(distnace2) < tol)
            cl_id = np.argmax(total_distance)

            if total_distance[cl_id] == 0:
                raise RuntimeError("The two points provided have to be on the same " + \
                                   " line (from inlet to outlet), and not at two different" + \
                                   " outlets")
            startID = ids1[cl_id]
            endID = ids2[cl_id]

    # Extract and spline a single line
    line_to_change = extract_single_line(centerline, cl_id, startID=startID, endID=endID)
    line_to_change = spline_centerline(line_to_change, nknots=25, isline=True)

    return line_to_change, region_points


def get_factor(line_to_change, method, beta, ratio, percentage, region_of_interest,
               region_points, stenosis_length, surface):
    """
    Compute the factor determining
    the change in radius, used to
    manipulate the Voronoi diagram.
    Only 80% of the geometry covered
    by the centerline is adjusted to
    leave the inlet, and the end of
    the area of interest, unchanged.

    Args:
        line_to_change (vtkPolyData): Centerline representing area of interest.
        beta (float): Factor deciding how area will change. Ignored if ratio is given.
        ratio (float): Desired ratio between min and max cross-sectional area.
        percentage (float): Desired increase/decrease in cross-sectional area.
        region_of_interest (str): Method for setting the region of interest.
        region_points (bool): List of points for the stenosis.
        stenosis_length (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Array to change the radius
    area = get_array("CenterlineSectionArea", line_to_change)

    # Safety smoothing, section area does not always work perfectly
    for i in range(2):
        area = gaussian_filter(area, 5)
    mean = np.mean(area)

    # Exclude first and last 5 % for some combinations of method an region_of_interest
    if region_of_interest in ["manuall", "commandline"] and method in ["area", "variation"]:
        k = int(round(factor_.shape[0] * 0.05, 0))
        l = area.shape[0] - k*2
    else:
        k = 0
        l = area.shape[0]
    trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(l).tolist() +
                        np.linspace(0, 1, k).tolist())

    # Get factor
    if method == "variation":
        if ratio is not None:
            # Inital guess
            R_old = area.max() / area.min()
            beta = 0.5 * math.log(ratio / R_old) / math.log(ratio) + 1

            # Parameters for algorithm
            R = 1e10
            a = 0
            b = 2
            sign = 0 if ratio > R_old else 1
            max_iter = 30
            iter = 0

            # Estimate beta with bisection method
            while abs(R - ratio) >= 0.001 and iter < max_iter:
                factor_ = (area / mean)**(beta-1)
                factor = factor_[:, 0]*(1-trans) + trans

                area_new = (np.sqrt(area[:, 0]/math.pi)*factor)**2 * math.pi
                R = area_new.max() / area_new.min()

                if R < ratio:
                    a = beta
                    beta = a + (b - a) / 2.
                else:
                    b = beta
                    beta = a + (b - a) / 2.

                iter += 1

            beta = beta - 1

        else:
            factor_ = ((area / mean)**beta)[:, 0]
            factor = factor_ * (1 - trans) + trans

    elif method == "stenosis":
        if len(region_points) == 3:
            t = np.linspace(0, np.pi, line_to_change.GetNumberOfPoints())
            factor = (1 - np.sin(t) * percentage * 0.01).tolist()
            factor = factor * (1 - trans) + trans

        elif len(region_points) == 6:
            length = get_curvilinear_coordinate(line_to_change)
            factor = area[0] + (area[-1] - area[0]) * (length / length.max())
            factor = factor * (1- trans) + trans

    # Increase or deacrease overall area by a percentage
    elif method == "area":
        factor_ = np.ones(len(trans)) * (1 + percentage*0.01)
        factor = factor_ * (1 - trans) + trans

    return factor


def change_area(voronoi, line_to_change, method, beta, ratio, percentage,
                region_of_interest, region_points, stenosis_length, surface):
    """
    Change the cross-sectional area of an input
    voronoi diagram along the corresponding area
    represented by a centerline.

    Args:
        voronoi (vtkPolyData): Voronoi diagram.
        tol (float): Tolerance factor.
        beta (float): Factor deciding how area will change. Ignored if ratio is given.
        ratio (float): Desired ratio between min and max cross-sectional area.
        percentage (float): Percentage the area of the geometry / stenosis is increase/decreased.
        region_of_interest (str): Method for setting the region of interest ['manuall' | 'commandline' | 'first_line']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        stenosis_length (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.

    Returns:
        newVoronoi (vtkPolyData): Manipulated Voronoi diagram.
    """
    arrayForTube = get_vtk_array("TubeRadius", 1, line_to_change.GetNumberOfPoints())

    # Note: If you are looking at the ICA or other vascular structure which can be rather
    # non-circular, please multiply MISR by a factor of, e.g. ~1.7. For arteries that are
    # very close, there is an oposite problem, and one can include points from the Voronoi
    # diagram belonging to other arteries. For robustness, the factor is now 1.05.
    MISR = get_array(radiusArrayName, line_to_change)*1.05
    for i in range(MISR.shape[0]):
        arrayForTube.SetTuple1(i, MISR[i])
    line_to_change.GetPointData().AddArray(arrayForTube)

    tubeFunction = vtkvmtk.vtkvmtkPolyBallLine()
    tubeFunction.SetInput(line_to_change)
    tubeFunction.SetPolyBallRadiusArrayName("TubeRadius")

    # Make a sphere at the end of the line
    pointID = line_to_change.GetNumberOfPoints()
    c1 = line_to_change.GetPoints().GetPoint(pointID - 1)
    c2 = line_to_change.GetPoints().GetPoint(pointID - 2)
    r = get_array(radiusArrayName, line_to_change)[-1]
    t = [c1[i] - c2[i] for i in range(len(c1))]
    lastSphere = vtk.vtkSphere()
    lastSphere.SetRadius(r * 1.5)
    lastSphere.SetCenter(c1)

    # Get factor
    factor = get_factor(line_to_change, method, beta, ratio, percentage,
                        region_of_interest, region_points, stenosis_length, surface)

    # Locator to find closest point on centerline
    locator = get_locator(line_to_change)
    M = line_to_change.GetNumberOfPoints()

    # Voronoi diagram
    N = voronoi.GetNumberOfPoints()
    newVoronoi = vtk.vtkPolyData()
    voronoiPoints = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, N)

    # Iterate through Voronoi diagram and manipulate
    # If inside MISR tube and inside plane, change r and move point.
    point = [0., 0., 0.]
    for i in range(N):
        voronoi.GetPoint(i, point)
        pointRadius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

        tubeValue = tubeFunction.EvaluateFunction(point)
        sphereValue = lastSphere.EvaluateFunction(point)
        voronoiVector = [point[j] - c1[j] for j in range(3)]
        vectorValue = vtk.vtkMath.Dot(voronoiVector, t)

        if (tubeValue <= 0.0) and not ((sphereValue < 0.0) and (vectorValue < 0.0)):
            tmp_ID = locator.FindClosestPoint(point)
            v1 = np.asarray(line_to_change.GetPoint(tmp_ID)) - np.asarray(point)
            v2 = v1 * (1 - factor[tmp_ID])
            point = (np.asarray(point) + v2).tolist()

            # Change radius
            pointRadius = pointRadius*factor[tmp_ID]

        voronoiPoints.InsertNextPoint(point)
        radiusArray.SetTuple1(i, pointRadius)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i)

    newVoronoi.SetPoints(voronoiPoints)
    newVoronoi.SetVerts(cellArray)
    newVoronoi.GetPointData().AddArray(radiusArray)

    return newVoronoi


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
    required = parser.add_argument_group('required named arguments')

    # Required arguments
    required.add_argument('-i', '--ifile', type=str, default=None, required=True,
                          help="Path to the surface model")
    # Output filename
    required.add_argument("-o", "--ofile", type=str, default=None, required=True,
                          help="Relative path to the output surface. The default folder is" + \
                               " the same as the input file, and a name with a combination of the" + \
                               " parameters.")


    # General arguments
    parser.add_argument("-m", "--method", type=str, default="variation",
                        choices=["variation", "stenosis", "area"],
                        help="Methods for manipulating the area in the region of interest:" + \
                             "\n1) 'variation' will increase or decrease the changes in area" + \
                             " along the centerline of the region of interest." + \
                             "\n2) 'stenosis' will create or remove a local narrowing of the" + \
                             " surface. If two points is provided, the area between these" + \
                             " two points will be linearly interpolated to remove the narrowing." + \
                             " If only one point is provided it is assumed to be the center of" + \
                             " the stenosis. The new stenosis will have a sin shape, however, any" + \
                             " other shape may be easly implemented." + \
                             "\n3) 'area' will inflate or deflate the area in the region of" + \
                             " interest.")
    parser.add_argument('-s', '--smooth', type=str2bool, default=True,
                        help="Smooth the voronoi diagram, default is False")
    parser.add_argument('-f', '--smooth_factor', type=float, default=0.25,
                        help="If smooth option is true then each voronoi point" + \
                             " that has a radius less then MISR*(1-smooth_factor) at" + \
                             " the closest centerline point is removed.")
    parser.add_argument("-n", "--no_smooth", type=bool, default=False,
                        help="If true and smooth is true the user, if no_smooth_point is" + \
                             " not given, the user can provide points where the surface not will" + \
                             " be smoothed.")
    parser.add_argument("--no_smooth_point", nargs="+", type=float, default=None,
                        help="If model is smoothed the user can manually select points on" + \
                             " the surface that will not be smoothed. A centerline will be" + \
                             " created to the extra point, and the section were the centerline" + \
                             " differ from the other centerlines will be keept un-smoothed. This" + \
                             " can be practicle for instance when manipulating geometries" + \
                             " with aneurysms")
    parser.add_argument("-b", "--poly-ball-size", nargs=3, type=int, default=[120, 120, 120],
                        help="The size of the poly balls that will envelope the new" + \
                             " surface. The default value is 120, 120, 120. If two tubular" + \
                             " structures are very close compared to the bounds, the poly ball" + \
                             " size should be adjusted. For quick proto typing we" + \
                             " recommend ~100 in all directions, but >250 for a final " + \
                             " surface.", metavar="size")
    parser.add_argument("-r", "--region-of-interest", type=str, default="manuall",
                        choices=["manuall", "commandline", "first_line"],
                        help="The method for defining the region to be changed. There are" + \
                             " three options: 'manuall', 'commandline', 'first_line'. In" + \
                             " 'manuall' the user will be provided with a visualization of the" + \
                             " input surface, and asked to provide an end and start point of the" + \
                             " region of interest. Note that not all algorithms are robust over" + \
                             " bifurcations. If 'commandline' is provided, then '--region-points'" + \
                             " is expected to be provided. Finally, if 'first_line' is given, the" + \
                             " line from the inlet (largest opening) to the first bifurcation" + \
                             " will be altered, not that method='stenosis' can not be used" + \
                             " with 'first_line'.")
    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" + \
                             " argument have to be given. The method expects two points" + \
                             " which defines the start and end of the region of interest. If" + \
                             " 'method' is set to stenosis, then one point can be provided as well," + \
                             " which is assumbed to be the center of a new stenosis." + \
                             " Example providing the points (1, 5, 1) and (2, 4, 3):" + \
                             " --stenosis-points 1 5 1 2 4 3")

    # "Variation" argments
    parser.add_argument('--beta', type=float, default=0.5,
                        help="For method=variation: The new voronoi diagram is computed as" + \
                             " (A/A_mean)**beta*r_old, over the respective area. If beta <" + \
                             " -1 the geometry will have less area variation, and" + \
                             " if beta > 1, the variations in area will increase")
    parser.add_argument("--ratio", type=float, default=None,
                        help="For method=variation: " + \
                             " Target ratio of A_max/A_min, when this is given beta will be ignored" + \
                             " and instead approximated to obtain the target ratio")

    # "Stenosis" argument
    parser.add_argument("--stenosis-length", type=float, default=2.0, metavar="length",
                        help="For method=stenosis: The length of the area " + \
                             " affected by a stenosis relative to the minimal inscribed" + \
                             " sphere radius of the selected point. Default is 2.0.")

    # "area" / "stenosis" argument
    parser.add_argument("--percentage", type=float, default=50.0,
                        help="Percentage the" + \
                             " area of the geometry is increase/decreased overall or only" + \
                             " stenosis")

    # Outputfile argument
    args = parser.parse_args()

    if args.method == "stenosis" and args.region_of_interest == "first_line":
        raise ValueError("Can not set region of interest to 'first_line' when creating or" + \
                         " removing a stenosis")

    if args.method == "variation" and args.ratio is not None and args.beta != 0.5:
        print("WARNING: The beta value you provided will be ignored, using ration instead.")

    if args.region_points is not None:
        if len(args.region_points) % 3 != 0 or len(args.region_points) > 6:
            raise ValueError("ERROR: Please provide region point(s) as a multiple of 3, and maximum" + \
                             " two points.")

    if args.no_smooth_point is not None and len(args.no_smooth_point):
        if len(args.no_smooth_point) % 3 != 0:
            raise ValueError("ERROR: Please provide the no smooth point(s) as a multiple" + \
                             " of 3.")


    return dict(input_filepath=args.ifile, method=args.method, smooth=args.smooth,
                smooth_factor=args.smooth_factor, beta=args.beta,
                region_of_interest=args.region_of_interest,
                region_points=args.region_points, ratio=args.ratio,
                stenosis_length=args.stenosis_length,
                percentage=args.percentage, output_filepath=args.ofile,
                poly_ball_size=args.poly_ball_size, no_smooth=args.no_smooth,
                no_smooth_point=args.no_smooth_point)

if __name__ == '__main__':
    area_variations(**read_command_line())

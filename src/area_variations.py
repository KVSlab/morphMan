from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelextrema
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from os import path, listdir

# Local import
from common import *

def area_variations(input_filepath, method, smooth, smooth_factor, no_smooth,
                    no_smooth_point, beta, ratio, stenosis_size, stenosis_points,
                    percentage, output_filepath, poly_ball_size):
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
        beta (float): Factor determining how much the geometry will differ.
        ratio (float): Target ratio, A_max / A_min. Beta is ignored (and estimated) if given.
        percentage (float): Percentage the area of the geometry / stenosis is increase/decreased.
        stenosis_size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        stenosis_points (float): Points of the stenosis. If one point is provided it is
        assumed to be the center of the new stenosis, and if two points are provided then
        they are assumbed to be the start and endpoint of an existing stenosis.

    Returns:
        length (ndarray): Array of abscissa coordinates.
    Returns:
        area (ndarray): Array of cross-section area along the abscissa.
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

    # Output path of manipulated surface
    if output_filepath is None:
        s = "_" + method
        s += "" if not smooth else "_smooth"
        if method == "variation":
            s += "_beta_" + str(beta) if beta is not None else "_ratio_" + str(ratio)
        elif method == "area":
            s += "_area_" + str(percentage)
        elif method == "stenosis":
            if stenosis_points is None:
                s += "_manual_size_{}_precentege_{}".format(stenosis_size, percentage)
            if stenosis_points is not None and len(stenosis_points) == 3:
                s += "_length_{}_narrowing_{}".format(stenosis_size, percentage)
            elif stenosis_points is not None and len(stenosis_points) == 6:
                s += "_remove_stenosis"
        new_surface_path =  path.join(base_path + "{}.vtp".format(s))
    else:
        new_surface_path = output_filepath

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Import centerline
    inlet, outlets = get_centers(surface, base_path)
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=0.1,
                                                         smooth=False, base_path=base_path)
    # Smooth voronoi diagram
    if smooth:
        voronoi = prepare_voronoi_diagram(surface, capped_surface, centerlines, base_path,
                                          smooth, smooth_factor, no_smooth,
                                          no_smooth_point, voronoi, pole_ids)
    tolerance = get_tolerance(centerlines)

    # Spline centerline and compute cross-sectional areas along line
    if not path.exists(centerline_area_spline_path):
        centerline_splined = get_line_to_change(centerlines, tolerance)
        write_polydata(centerline_splined, centerline_spline_path)

        centerline_area, centerline_area_sections = vmtk_compute_centerline_sections(surface,
                                                                              centerline_splined)
        write_polydata(centerline_area, centerline_area_spline_path)
        write_polydata(centerline_area_sections, centerline_area_spline_sections_path)
    else:
        centerline_area = read_polydata(centerline_area_spline_path)

    # Manipulate the voronoi diagram
    print("Change Voronoi diagram")
    newvoronoi = change_area(voronoi, centerline_area, method, beta, ratio, percentage,
                             stenosis_points, stenosis_size, surface)
    write_polydata(newvoronoi, voronoi_new_path)

    # Make new surface
    print("Create surface")
    new_surface = create_new_surface(newvoronoi, poly_ball_size=poly_ball_size)

    print("Cleaning surface for output")
    new_surface = prepare_surface_output(new_surface, surface, centerlines,
                                         new_surface_path, test_merge=True)

    print("Write surface to: {}".format(new_surface_path))
    write_polydata(new_surface, new_surface_path)


def get_line_to_change(centerline, tol):
    """
    Extract and spline part of centerline
    within the geometry where
    area variations will be performed.

    Args:
        centerline (vtkPolyData): Centerline in geometry.
        tol (float): Tolerance parameter.

    Returns:
        line_to_change (vtkPolyData): Part of centerline.
    """
    line2 = extract_single_line(centerline, 0)
    numberOfPoints2 = line2.GetNumberOfPoints()

    n = 2
    pointIDs = []
    # Iterate through lines and find diverging point
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

    pointID = min(pointIDs)

    # Extract and spline a single line
    line_to_change = extract_single_line(centerline, 0, endID=pointID)
    line_to_change = spline_centerline(line_to_change, nknots=25, isline=True)

    return line_to_change


def get_factor(line_to_change, method, beta, ratio, percentage, stenosis_points, stenosis_size, surface):
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
        stenosis_points (bool): List of points for the stenosis.
        stenosis_size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Array to change the radius
    area = get_array("CenterlineSectionArea", line_to_change)
    for i in range(2):
        area = gaussian_filter(area, 5)
    mean = np.mean(area)

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

            # Exclude first and last 10 %
            area_ = area[int(area.shape[0]*0.02):-int(area.shape[0]*0.02)]

            # Estimate beta and find factor
            while abs(R - ratio) >= 0.001 and iter < max_iter:
                factor_ = (area / mean)**(beta-1)
                k = int(round(factor_.shape[0] * 0.10, 0))
                l = factor_.shape[0] - k*2
                trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(l).tolist() +
                                   np.linspace(0, 1, k).tolist())
                factor = factor_[:, 0]*(1-trans) + trans

                area_new = (np.sqrt(area[:, 0]/math.pi)*factor)**2 * math.pi
                R = area_new.max() / area_new.min()

                print(("R now: {:4f}  Want: {:2f}  R old: {:6f}".format(R, ratio, R_old)))
                if R < ratio:
                    a = beta
                    beta = a + (b - a) / 2.
                else:
                    b = beta
                    beta = a + (b - a) / 2.

                iter += 1

            beta = beta - 1

        else:
            factor_ = (area / mean)**beta

            # A linear transition of the old and new geometry
            k = int(round(factor_.shape[0] * 0.10, 0))
            l = factor_.shape[0] - k*2
            trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(l).tolist() +
                               np.linspace(0, 1, k).tolist())
            factor = factor_[:, 0]*(1-trans) + trans


    elif method == "stenosis":
        if stenosis_points is not None:
            stenosisPointId = vtk.vtkIdList()
            first = True
            while stenosisPointId.GetNumberOfIds() not in [1, 2]:
                if not first:
                    print("Please provide only one or two points, try again.")
                # Select point on surface
                seed_selector = vmtkPickPointSeedSelector()
                seed_selector.SetSurface(surface)
                seed_selector.Mode = "stenosis"
                seed_selector.Execute()
                stenosisPointId = seed_selector.GetTargetSeedIds()
                first = False

            if stenosisPointId.GetNumberOfIds() == 1:
                print("One point detected. Creating stenosis.")
                stenosisPoint = surface.GetPoint(stenosisPointId.GetId(0))
                factor = create_stenosis(line_to_change, stenosisPoint, area, stenosis_size, percentage)

            elif stenosisPointId.GetNumberOfIds() == 2:
                print("Two points detected. Removeing stenosis between the two points.")
                stenosisPoint1 = surface.GetPoint(stenosisPointId.GetId(0))
                stenosisPoint2 = surface.GetPoint(stenosisPointId.GetId(1))
                factor = remove_stenosis(line_to_change, stenosisPoint1, stenosisPoint2, area)

        else:
            if len(stenosis_points) == 6:
                factor = remove_stenosis(line_to_change, stenosis_points[:3],
                                         stenosis_points[3:], area)
            else:
                factor = create_stenosis(line_to_change, stenosis_points, area, stenosis_size, percentage)

    elif method == "area":
        # Increase or deacrease overall area by a percentage
        k = int(round(area.shape[0] * 0.10, 0))
        l = area.shape[0] - k*2
        trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(l).tolist() +
                           np.linspace(0, 1, k).tolist())
        factor_ = np.ones(len(trans)) * (1 + percentage*0.01)
        factor = factor_*(1-trans) + trans

    return factor


def create_stenosis(line_to_change, stenosisPoint, area, stenosis_size, percentage):
    """
    Creates a stenosis along the vessel based
    on a user selected point representing
    the center of the stenosis.

    Args:
        line_to_change (vtkPolyData): Centerline representing area of interest.
        stenosisPoint (tuple): Center of stenosis.
        area (ndarray): Array of cross-section area along the abscissa.
        stenosis_size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        percentage (float): Desired increase/decrease in cross-sectional area.

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Find closest point along centerline
    locator = get_locator(line_to_change)
    stenosisLinePointId = int(locator.FindClosestPoint(stenosisPoint))
    stenosisLinePoint = line_to_change.GetPoint(stenosisLinePointId)

    # Find start and stop based on MISR
    stenosisLinePointMISR = get_array(radiusArrayName, line_to_change)[stenosisLinePointId]
    tolerance = stenosisLinePointMISR * stenosis_size
    startId = 0
    stopId = line_to_change.GetNumberOfPoints()-1
    for i in range(stenosisLinePointId + 1, 0, -1):
        stenosisPoint_tmp = line_to_change.GetPoint(i)
        if distance(np.asarray(stenosisPoint_tmp), np.asarray(stenosisLinePoint)) > tolerance:
            startId = i
            break

    for i in range(stenosisLinePointId + 1, line_to_change.GetNumberOfPoints(), 1):
        stenosisPoint_tmp = line_to_change.GetPoint(i)
        if distance(np.asarray(stenosisPoint_tmp), np.asarray(stenosisLinePoint)) > tolerance:
            stopId = i
            break

    # Subdivide centerline into sections
    trans = subdivide_centerline(area.shape[0], startId, stopId, segment_length=0.9)

    # Define sine profile
    t = np.linspace(0, np.pi, stopId-startId)
    sine_profile = (np.sin(t)*percentage*0.01 + 1.0).tolist()
    factor_ = np.asarray(np.ones(startId).tolist() + sine_profile + np.ones(area.shape[0]-stopId).tolist())
    factor = factor_*(1-trans) + trans

    return factor


def subdivide_centerline(size, startId, stopId, segment_length=1.0):
    """
    Create a linear transition between original
    and manipulated geometry, by dividing
    centerline into sections.

    Args:
        size (int): Number of points along centerline.
        startId (int): ID at starting point.
        stopId (int): ID at stopping point.
        segment_length (float): Fraction of centerline used.

    Returns:
        trans (ndarray): Linear transition between original and new geomtry.
    """
    l = stopId - startId
    diff = size - stopId
    end_l = int(diff * segment_length)
    endmid_l = diff - end_l
    start_l = int(startId * segment_length)
    startmid_l = startId - start_l

    start = np.asarray(np.ones(start_l)).tolist()
    startmid = np.asarray(np.linspace(1, 0, startmid_l)).tolist()
    mid = np.zeros(l).tolist()
    endmid = np.linspace(0, 1, endmid_l).tolist()
    end = np.asarray(np.ones(end_l)).tolist()
    trans = np.asarray(start + startmid + mid + endmid + end)

    return trans


def remove_stenosis(line_to_change, point1, point2, area):
    """
    Removes a stenosis along the vessel based
    on a two user selected point representing
    the boundary of the stenosis.

    Args:
        line_to_change (vtkPolyData): Centerline representing area of interest.
        point 1 (tuple): Start point of stenosis.
        point 2 (tuple): End point of stenosis.
        area (ndarray): Array of cross-section area along the abscissa.

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Find closest point along centerline and sort points
    locator = get_locator(line_to_change)
    lineId1 = int(locator.FindClosestPoint(point1))
    lineId2 = int(locator.FindClosestPoint(point2))
    if lineId1 < lineId2:
        startId = lineId1
        stopId = lineId2
    else:
        startId = lineId2
        stopId = lineId1

    # Get areas at start and end
    startPoint = line_to_change.GetPoint(startId)
    endPoint = line_to_change.GetPoint(stopId)
    startMISR = get_array(radiusArrayName, line_to_change)[startId]
    endMISR = get_array(radiusArrayName, line_to_change)[stopId]

    def new_area(s):
        tmpId = int(s)
        tmprad = get_array(radiusArrayName, line_to_change)[tmpId]
        s = (s - startId) / (stopId - startId)
        return (startMISR * (1 - s) + s * endMISR) / tmprad

    # Subdivide centerline into sections
    trans = subdivide_centerline(area.shape[0], startId, stopId, segment_length=1.0)

    # Define stenosis area and factor
    t = np.linspace(startId, stopId, stopId-startId)
    fixed_area = []
    for t_ in t:
        fixed_area.append(new_area(t_))
    factor_ = np.asarray(np.ones(startId).tolist() + fixed_area + np.ones(area.shape[0]-stopId).tolist())
    factor = factor_*(1-trans) + trans

    return factor


def change_area(voronoi, line_to_change, method, beta, ratio, percentage, stenosis_points,
                stenosis_size, surface):
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
        stenosis_points (list): List of points for the stenosis.
        stenosis_size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.

    Returns:
        newVoronoi (vtkPolyData): Manipulated Voronoi diagram.
    """
    arrayForTube = get_vtk_array("TubeRadius", 1, line_to_change.GetNumberOfPoints())
    MISR = get_array(radiusArrayName, line_to_change)*1.7
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
    factor = get_factor(line_to_change, method, beta, ratio, percentage, stenosis_points,
                        stenosis_size, surface)

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

    # General arguments
    parser.add_argument("-m", "--method", type=str, default="variation",
                        choices=["variation", "stenosis", "area"], metavar="method",
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
    parser.add_argument('-s', '--smooth', type=bool, default=False,
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

    # "Variation" argments
    parser.add_argument('--beta', type=float, default=0.5,
                        help="For method=variation: The new voronoi diagram is computed as" + \
                             " (A/A_mean)**beta*r_old, over the respective area. If beta <" + \
                             " -1 the geometry will have less area variation, and" + \
                             " if beta > 1, the variations in area will increase")
    parser.add_argument("--ratio", type=float, default=None, help="For method=variation: " + \
                        "Target ratio of A_max/A_min, when this is given beta will be ignored" + \
                       " and instead approximated to obtain the target ratio")

    # "Stenosis" argument
    parser.add_argument("--stenosis-size", type=float, default=2.0, metavar="size",
                        help="For method=stenosis: The length of the area " + \
                        " affected by a stenosis. Default is 2.0 times the minimal" + \
                        " inscribed spehere radius of the selected point.")
    parser.add_argument("--stenosis-points", nargs="+", type=float, default=None, metavar="points",
                        help="If the points are not provided, the user will be asked to" + \
                        " give them provide the points by clicking on the input surface." + \
                        " If one point is provided it is assumed to be the center of the" + \
                        " stenosis, and if two points are provided it is assumed to be" + \
                        " the start and end of an existing stenosis. Example providing two" + \
                        " points: --stenosis-points 1.0 1.0 1.0 2.0 2.0 2.0")

    # "area" / "stenosis" argument
    parser.add_argument("--percentage", type=float, default=50.0, help="Percentage the" + \
                        " area of the geometry is increase/decreased overall or only stenosis")

    # Output filename
    parser.add_argument("-o", "--ofile", type=str, default=None,
                        help="Relative path to the output surface. The default folder is" + \
                        " the same as the input file, and a name with a combination of the" + \
                        " parameters.")

    # Outputfile argument
    args = parser.parse_args()

    if args.method == "variation" and args.ratio is not None and args.beta != 0.5:
        print("WARNING: The beta value you provided will be ignored, using ration instead.")

    if args.stenosis_points is not None and len(args.stenosis_points):
        if len(args.stenosis_points) % 3 != 0 and len(args.stenosis_points) < 7:
            raise ValueError("ERROR: Please provide a stenosis point as a multiple of 3, and maximum" + \
                             " two points")

    if args.no_smooth_point is not None and len(args.no_smooth_point):
        if len(args.no_smooth_point) % 3 != 0:
            raise ValueError("ERROR: Please provide the no smooth point(s) as a multiple of 3")


    return dict(input_filepath=args.ifile, method=args.method, smooth=args.smooth,
                smooth_factor=args.smooth_factor, beta=args.beta,
                ratio=args.ratio, stenosis_size=args.stenosis_size, stenosis_points=args.stenosis_points,
                percentage=args.percentage, output_filepath=args.ofile,
                poly_ball_size=args.poly_ball_size, no_smooth=args.no_smooth,
                no_smooth_point=args.no_smooth_point)

if __name__ == '__main__':
    area_variations(**read_command_line())

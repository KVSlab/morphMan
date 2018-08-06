from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelextrema
from argparse import ArgumentParser
from os import path, listdir

# Local import
from common import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *


def read_command_line():
    """
    Read arguments from commandline
    """
    parser = ArgumentParser()

    parser.add_argument('-c', '--case', type=str, default=None, help="Case")
    parser.add_argument('-d', '--dir_path', type=str, default=".", help="Path")
    parser.add_argument('--s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram", metavar="smooth")
    parser.add_argument('--beta', type=float, default=0.5,
            help="The new voronoi diagram is computed as (A/mean)**beta*r_old," + \
            " over the respective area. If beta < -1 the geometry will be more even, and" + \
            " if beta > 1, the differences in the geometry will be larger")
    parser.add_argument("--ratio", type=float, default=None, help="Wanted ratio" + \
                       " A_max/A_min, when this is given beta will be ignored" + \
                       " and beta computed such that this will (approxematly)" + \
                       " be the result")
    parser.add_argument('-sf', '--smooth_factor', type=float, default=0.25,
                         help="If smooth option is true then each voronoi point" + \
                         " that has a radius less then MISR*(1-smooth_factor) at" + \
                         " the closest centerline point is removes", metavar="smoothening_factor")
    parser.add_argument("--percentage", type=float, default=None, help="Percentage the" + \
                        " area of the geometry is increase/decreased overall or only stenosis")
    parser.add_argument("--stenosis", type=bool, default=False, help="Creates a user selected" + \
                        " stenosis in in the geometry. Area is decreased by the parameter 'precentage'.")
    parser.add_argument("--size", type=float, default=2.0, help="The length of the area " + \
                        " affected by a stenosis. Default is MISR x 2.0 of selected point.")
    parser.add_argument("--noise", type=bool, default=False, help="Adds noise to surface area.")
    parser.add_argument("--stats", type=bool, default=False,
                        help="Collect stats")

    args = parser.parse_args()

    if args.ratio is not None and args.beta != 0.5:
        print("WARNING: The beta value you provided will be ignored.")

    return args.s, args.beta, args.stats, args.dir_path, args.case, args.ratio, \
           args.percentage, args.stenosis, args.size, args.noise, args.smooth_factor


def get_stats(centerline_area, folder, centerline):
    """
    Compute and write all statistical
    paramaters of the manipulation.

    Args:
        centerline_area (vtkPolyData): Centerline object including area.
        folder (str): Path to case folder.
        centerline (vtkPolyData): Centerline object.

    Returns:
        length (ndarray): Array of abscissa coordinates.
    Returns:
        area (ndarray): Array of cross-section area along the abscissa.
    """

    area = get_array("CenterlineSectionArea", centerline_area)
    MISR_ = get_array(radiusArrayName, centerline)**2*math.pi
    MISR = np.array([MISR_[i] for i in range(MISR_.shape[0] - 1, -1, -1)])[:area.shape[0]]
    length = get_curvilinear_coordinate(centerline_area)

    for i in range(2):
        area = gaussian_filter(area, 5)
        MISR = gaussian_filter(MISR, 5)

    # Area shape - "circleness"
    circleness = MISR / area
    max_circleness = circleness.max()
    min_circleness = circleness.min()
    mean_circleness = circleness.mean()

    # Local extremas, ignore min or max on boundary
    local_min_MISR_ID = argrelextrema(MISR, np.less)[0]
    local_max_MISR_ID = argrelextrema(MISR, np.greater)[0]
    local_min_area_ID = argrelextrema(area, np.less)[0]
    local_max_area_ID = argrelextrema(area, np.greater)[0]

    local_min_MISR = MISR[local_min_MISR_ID]
    local_max_MISR = MISR[local_max_MISR_ID]
    local_min_area = area[local_min_area_ID]
    local_max_area = area[local_max_area_ID]

    global_min_MISR = local_min_MISR.min()
    global_max_MISR = local_max_MISR.max()
    global_min_area = local_min_area.min()
    global_max_area = local_max_area.max()

    mean_area = area.mean()
    mean_MISR = MISR.mean()

    number_of_max = local_max_area.shape[0]

    # Min max derived parameters
    length_min_max = abs(length[(area == global_min_area).nonzero()[0]] - \
                         length[(area == global_max_area).nonzero()[0]])[0]

    max_mean_ratio_area = global_max_area / mean_area
    min_mean_ratio_area = global_min_area / mean_area
    max_mean_ratio_MSIR = global_max_MISR / mean_MISR
    min_mean_ratio_MISR = global_min_MISR / mean_MISR

    # Global and local disent
    global_min_max_disent = abs(math.sqrt(global_max_area)/math.pi -
                                math.sqrt(global_min_area) / math.pi) / \
                            length_min_max
    local_max_stepest_disent = 0

    if length[(area == local_min_area[0]).nonzero()[0]] > \
       length[(area == local_max_area[0]).nonzero()[0]]:
        start = 1
    else:
        start = 0

    N = min(number_of_max, local_min_area.shape[0] - start)
    for i in range(N):
        min_ = local_min_area[start + i]
        max_ = local_max_area[i]

        h = math.sqrt(max_)/math.pi - math.sqrt(min_) / math.pi
        l = abs(length[(area == max_).nonzero()[0]] - length[(area == min_).nonzero()[0]])
        if h/l > local_max_stepest_disent:
            local_max_stepest_disent = h / l

    # Max point disent (max |derivative|)
    knots = 40
    t = np.linspace(length[0], length[-1], knots+2)[1:-1]
    spline_area = splrep(length, area, k=4, t=t)
    spline_area_ = splev(length, spline_area)
    darea_dx = splev(length, spline_area, der=1)
    max_derivative = abs(darea_dx).max()

    stats = {"max_derivative": max_derivative,
             "local_max_stepest_disent": local_max_stepest_disent,
             "max_mean_ratio_area": max_mean_ratio_area,
             "min_mean_ratio_area": min_mean_ratio_area,
             "mean_area": mean_area,
             "max_min_ratio_area": global_max_area / global_min_area,
             "length_min_max": length_min_max,
             "global_min_area": global_min_area,
             "global_max_area": global_max_area,
             "max_circleness": max_circleness,
             "min_circleness": min_circleness,
             "mean_circleness": mean_circleness,
             "number_of_max": number_of_max}

    write_parameters(stats, folder)

    return length, area


def get_line_to_change(centerline, tol):
    """
    Extract and spline part of centerline
    within the geometry where
    area variations will be performed.

    Args:
        centerline (vtkPolyData): Centerline in geometry.
        tol (float): Tolerance parameter.

    Returns:
        lineToChange (vtkPolyData): Part of centerline.
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
    lineToChange = extract_single_line(centerline, 0, endID=pointID)
    lineToChange = spline_centerline(lineToChange, nknots=25, isline=True, get_stats=False)

    return lineToChange


def get_factor(lineToChange, beta, ratio, percentage, stenosis, stenosis_size, add_noise, surface):
    """
    Compute the factor determining
    the change in radius, used to
    manipulate the Voronoi diagram.
    Only 80% of the geometry covered
    by the centerline is adjusted to
    leave the inlet, and the end of
    the area of interest, unchanged.

    Args:
        lineToChange (vtkPolyData): Centerline representing area of interest.
        beta (float): Factor deciding how area will change. Ignored if ratio is given.
        ratio (float): Desired ratio between min and max cross-sectional area.
        percentage (float): Desired increase/decrease in cross-sectional area.
        stenosis (bool): Creates or removes a user selected stenosis in the geometry if True.
        stenosis_size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        add_noise (bool): Adds noise to surface if true.

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Array to change the radius
    area = get_array("CenterlineSectionArea", lineToChange)
    for i in range(2):
        area = gaussian_filter(area, 5)
    mean = np.mean(area)

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

    elif stenosis is True:
        # Select point on surface
        SeedSelector = vmtkPickPointSeedSelector()
        SeedSelector.SetSurface(surface)
        SeedSelector.Mode = "stenosis"
        SeedSelector.Execute()
        stenosisPointId = SeedSelector.GetTargetSeedIds()
        if stenosisPointId.GetNumberOfIds() == 1:
            if percentage is None:
                percentage = -40
                print("No increase/decrease factor given. Reducing max area by 40%")
            print("One point detected. Creating stenosis.")
            stenosisPoint = surface.GetPoint(stenosisPointId.GetId(0))
            factor = create_stenosis(lineToChange, stenosisPoint, area, stenosis_size, percentage)

        elif stenosisPointId.GetNumberOfIds() == 2:
            print("Two points detected. Removes stenosis between the two points.")
            stenosisPoint1 = surface.GetPoint(stenosisPointId.GetId(0))
            stenosisPoint2 = surface.GetPoint(stenosisPointId.GetId(1))
            factor = remove_stenosis(lineToChange, stenosisPoint1, stenosisPoint2, area)

    elif percentage is not None:
        # Increase or deacrease overall area by a percentage
        k = int(round(area.shape[0] * 0.10, 0))
        l = area.shape[0] - k*2
        trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(l).tolist() + np.linspace(0, 1, k).tolist())
        factor_ = np.ones(len(trans)) * (1 + percentage*0.01)
        factor = factor_*(1-trans) + trans

    else:
        factor_ = (area / mean)**beta
        # A linear transition of the old and new geometry
        k = int(round(factor_.shape[0] * 0.10, 0))
        l = factor_.shape[0] - k*2
        trans = np.asarray(np.linspace(1, 0, k).tolist() + np.zeros(l).tolist() +
                           np.linspace(0, 1, k).tolist())
        factor = factor_[:, 0]*(1-trans) + trans

    return factor


def create_stenosis(lineToChange, stenosisPoint, area, stenosis_size, percentage):
    """
    Creates a stenosis along the vessel based
    on a user selected point representing
    the center of the stenosis.

    Args:
        lineToChange (vtkPolyData): Centerline representing area of interest.
        stenosisPoint (tuple): Center of stenosis.
        area (ndarray): Array of cross-section area along the abscissa.
        stenosis_size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        percentage (float): Desired increase/decrease in cross-sectional area.

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Find closest point along centerline
    locator = get_locator(lineToChange)
    stenosisLinePointId = int(locator.FindClosestPoint(stenosisPoint))
    stenosisLinePoint = lineToChange.GetPoint(stenosisLinePointId)

    # Find start and stop based on MISR
    stenosisLinePointMISR = get_array(radiusArrayName, lineToChange)[stenosisLinePointId]
    tolerance = stenosisLinePointMISR * stenosis_size
    startId = 0
    stopId = lineToChange.GetNumberOfPoints()-1
    for i in range(stenosisLinePointId + 1, 0, -1):
        stenosisPoint_tmp = lineToChange.GetPoint(i)
        if distance(np.asarray(stenosisPoint_tmp), np.asarray(stenosisLinePoint)) > tolerance:
            startId = i
            break

    for i in range(stenosisLinePointId + 1, lineToChange.GetNumberOfPoints(), 1):
        stenosisPoint_tmp = lineToChange.GetPoint(i)
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


def remove_stenosis(lineToChange, point1, point2, area):
    """
    Removes a stenosis along the vessel based
    on a two user selected point representing
    the boundary of the stenosis.

    Args:
        lineToChange (vtkPolyData): Centerline representing area of interest.
        point 1 (tuple): Start point of stenosis.
        point 2 (tuple): End point of stenosis.
        area (ndarray): Array of cross-section area along the abscissa.

    Returns:
        factor (float): Factor determining the change in radius.
    """
    # Find closest point along centerline and sort points
    locator = get_locator(lineToChange)
    lineId1 = int(locator.FindClosestPoint(point1))
    lineId2 = int(locator.FindClosestPoint(point2))
    if lineId1 < lineId2:
        startId = lineId1
        stopId = lineId2
    else:
        startId = lineId2
        stopId = lineId1

    # Get areas at start and end
    startPoint = lineToChange.GetPoint(startId)
    endPoint = lineToChange.GetPoint(stopId)
    startMISR = get_array(radiusArrayName, lineToChange)[startId]
    endMISR = get_array(radiusArrayName, lineToChange)[stopId]

    def new_area(s):
        tmpId = int(s)
        tmprad = get_array(radiusArrayName, lineToChange)[tmpId]
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


def change_area(voronoi, lineToChange, beta, ratio, percentage, stenosis, stenosis_size, add_noise, surface):
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
        stenosis (bool): Creates or removes a user selected stenosis in the geometry if True.
        stenosis_size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        add_noise (bool): Adds noise to surface if true.

    Returns:
        newVoronoi (vtkPolyData): Manipulated Voronoi diagram.
    """
    arrayForTube = get_vtk_array("TubeRadius", 1, lineToChange.GetNumberOfPoints())
    MISR = get_array(radiusArrayName, lineToChange)*1.7
    for i in range(MISR.shape[0]):
        arrayForTube.SetTuple1(i, MISR[i])
    lineToChange.GetPointData().AddArray(arrayForTube)

    tubeFunction = vtkvmtk.vtkvmtkPolyBallLine()
    tubeFunction.SetInput(lineToChange)
    tubeFunction.SetPolyBallRadiusArrayName("TubeRadius")

    # Make a sphere at the end of the line
    pointID = lineToChange.GetNumberOfPoints()
    c1 = lineToChange.GetPoints().GetPoint(pointID - 1)
    c2 = lineToChange.GetPoints().GetPoint(pointID - 2)
    r = get_array(radiusArrayName, lineToChange)[-1]
    t = [c1[i] - c2[i] for i in range(len(c1))]
    lastSphere = vtk.vtkSphere()
    lastSphere.SetRadius(r * 1.5)
    lastSphere.SetCenter(c1)

    # Get factor
    factor = get_factor(lineToChange, beta, ratio, percentage, stenosis, stenosis_size, add_noise, surface)

    # Locator to find closest point on centerline
    locator = get_locator(lineToChange)
    M = lineToChange.GetNumberOfPoints()

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
            v1 = np.asarray(lineToChange.GetPoint(tmp_ID)) - np.asarray(point)

            # Add noise
            if 0.1*M < tmp_ID < 0.9*M and add_noise:
                v2 = v1 * (1 - factor[tmp_ID]*np.random.uniform(0.4, 1.2))
            else:
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


def area_variations(folder, beta, smooth, stats, r_change, percentage, stenosis,
                    stenosis_size, add_noise, smooth_factor):
    """
    Objective manipulation of area variation in
    patient-specific models of blood vessels.
    Manipulation is performed by looping over
    all Voronoi points relative to a selected centerline
    and change the corresponding radius.

    Args:
        folder (str): Path to directory with cases.
        beta (float): Factor determining how much the geometry will differ.
        smooth (bool): Determines Voronoi smoothing or not.
        stats (bool): Computes stats of area variations if True.
        r_change (float): Wanted ratio, A_max / A_min. Beta is ignored (and estimated) if given.
        percentage (float): Percentage the area of the geometry / stenosis is increase/decreased.
        stenosis (bool): Creates or removes a user selected stenosis in the geometry if True.
        stenosis_size (float): Length of affected stenosis area. Default is MISR x 2.0 of selected point.
        add_noise (bool): Adds noise to surface if true.
        smooth_factor (float): Smoothing factor used for voronoi diagram smoothing.

    Returns:
        length (ndarray): Array of abscissa coordinates.
    Returns:
        area (ndarray): Array of cross-section area along the abscissa.
    """
    # Input files
    model_path = path.join(folder, "surface", "model.vtp")
    centerlines_path = path.join(folder, "surface", "model_usr_centerline.vtp")
    voronoi_path = path.join(folder, "surface", "voronoi.vtp")

    # Output files
    voronoi_smoothed_path = path.join(folder, "surface", "voronoi_smoothed.vtp")
    voronoi_new_path = path.join(folder, "surface", "voronoi_area.vtp")
    centerline_area_path = path.join(folder, "surface", "centerline_area.vtp")
    centerline_area_spline_path = path.join(folder, "surface", "centerline_area_spline.vtp")
    centerline_area_spline_sections_path = path.join(folder, "surface", "centerline_area_sections.vtp")
    centerline_spline_path = path.join(folder, "surface", "centerline_spline.vtp")
    centerline_new_path = path.join(folder, "surface", "centerline_area_new.vtp")
    model_smoothed_path = path.join(folder, "surface", "model_smoothed.vtp")
    s = ""
    s += "" if not smooth else "_smooth"
    s += "" if r_change is not None or percentage is not None else "_%s" % beta
    s += "" if r_change is None else "_ratio%s" % r_change
    s += "" if not stenosis else "_stenosis_%smisr" % stenosis_size
    s += "" if not add_noise else "_noise"
    s += "" if percentage is None else "_%s" % percentage
    model_new_surface =  path.join(folder, "surface", "model_area%s.vtp" % s)
    model_new_surface_clean =  path.join(folder, "surface", "model_area%s_clean.vtp" % s)

    # Import centerline
    centerlines = make_centerline(model_path, centerlines_path, length=0.1, smooth=False)

    # Clean and capp / uncapp surface
    parameters = get_parameters(folder)
    surface, capped_surface = preare_surface(model_path, parameters)

    # Smooth voronoi diagram
    voronoi = prepare_voronoi_diagram(surface, model_smoothed_path, voronoi_path, voronoi_smoothed_path,
                                    smooth, smooth_factor, centerlines)

    # Tolerance for finding diverging point
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

    area = None
    length = None

    if stats:
        # Compute stats
        length, area = get_stats(centerline_area, folder, centerlines)
    else:
        # Change and compute the new voronoi diagram
        print("Change Voronoi diagram")
        newvoronoi = change_area(voronoi, centerline_area, beta, ratio, percentage,
                    stenosis, stenosis_size, add_noise, surface)

        print("Write Voronoi diagram")
        write_polydata(newvoronoi, voronoi_new_path)

        # Make new surface
        print("Create surface")
        new_surface = create_new_surface(newvoronoi)

        print("Write surface to: {}".format(model_new_surface.split("/")[-1]))
        # TODO: Add Automated clipping of newmodel
        new_surface = vmtk_surface_smoother(new_surface, method="laplace", iterations=100)
        new_surface = clean_and_check_surface(new_surface, centerlines_in_order,
                                        model_new_surface_clean, centerline_new_path)
        write_polydata(new_surface, model_new_surface)

    return length, area


if __name__ == '__main__':
    smooth, beta, stats, basefolder, case, ratio, percentage, stenosis, stenosis_size, add_noise, smooth_factor = read_command_line()
    folders = listdir(basefolder) if case is None else [case]
    for i, folder in enumerate(folders):
        if folder.startswith("P0"):
            print("Working on: {}".format(folder))
            case = path.join(basefolder, folder)
            area_variations(case, beta, smooth, stats, ratio, percentage, stenosis, stenosis_size, add_noise, smooth_factor)

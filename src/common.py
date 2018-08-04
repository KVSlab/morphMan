import vtk
import numpy as np
import numpy.linalg as la
import sys
import math
from vmtk import vtkvmtk, vmtkscripts
from os import path, listdir
from scipy.signal import argrelextrema, gaussian, resample
from scipy.interpolate import splrep, splev
from vmtkpointselector import *

# Global array names
radiusArrayName = 'MaximumInscribedSphereRadius'
parallelTransportNormalsArrayName = 'ParallelTransportNormals'
groupIDsArrayName = "GroupIds"
abscissasArrayName = 'Abscissas'
clippingArrayName = 'ClippingArray'
branchClippingArrayName = 'BranchClippingArray'
distanceToTubeArrayName = 'DistanceToTubeFunction'
closedArrayName = 'ClosedSection'
eikonalSolutionArrayName = 'EikonalSolutionArray'
edgeArrayName = 'EdgeArray'
edgePCoordArrayName = 'EdgePCoordArray'
costFunctionArrayName = 'CostFunctionArray'

# Options not available from commandline
divergingRatioToSpacingTolerance = 2.0
interpolationHalfSize = 3
voronoiCoreCutOffThreshold = 0.75
numberOfSplineAnalyzedPoints = 40
phiValues = [float(i) for i in range(2, 43, 2)]
thetaStep = 2.0


def read_polydata(filename):
    """
    Load the given file, and return a vtkPolyData object for it.

    Args:
        filename (str): Path to input file.

    Returns:
        polyData (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear∕vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Output data.
    """

    # Check if file exists
    if not path.exists(filename):
        raise RuntimeError("Could not find file: %s" % filename)

    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get reader
    if fileType == 'stl':
        reader = vtk.vtkSTLReader()
        reader.MergingOn()
    elif fileType == 'vtk':
        reader = vtk.vtkPolyDataReader()
    elif fileType == 'vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif fileType == 'vts':
        reader = vtk.vtkXMLStructuredGridReader()
    elif fileType == 'vtr':
        reader = vtk.vtkXMLRectilinearGridReader()
    elif fileType == 'vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif fileType == "vti":
        reader = vtk.vtkXMLImageDataReader()
    elif fileType == "tec":
        polyData = ReadTecplotSurfaceFile(filename)
        return polyData
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Read
    reader.SetFileName(filename)
    reader.Update()
    polyData = reader.GetOutput()

    return polyData


def write_polydata(input_data, filename):
    """
    Write the given input data based on the file name extension.

    Args:
        input_data (vtkSTL∕vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear∕vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Input data.
        filename (str): Save path location.
    """
    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get writer
    if fileType == 'stl':
        writer = vtk.vtkSTLWriter()
    elif fileType == 'vtk':
        writer = vtk.vtkPolyDataWriter()
    elif fileType == 'vts':
        writer = vtk.vtkXMLStructuredGridWriter()
    elif fileType == 'vtr':
        writer = vtk.vtkXMLRectilinearGridWriter()
    elif fileType == 'vtp':
        writer = vtk.vtkXMLPolyDataWriter()
    elif fileType == 'vtu':
        writer = vtk.vtkXMLUnstructuredGridWriter()
    elif fileType == "vti":
        writer = vtk.vtkXMLImageDataWriter()
    elif fileType == "tec":
        WriteTecplotSurfaceFile(input_data, filename)
        return
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Set filename and input
    writer.SetFileName(filename)
    writer.SetInputData(input_data)
    writer.Update()

    # Write
    writer.Write()


def make_voronoi_diagram(surface, filename):
    """
    Creates a surface model's
    coresponding voronoi diagram

    Args:
        surface (vtkPolyData): Surface model
        filename (str): Path where voronoi diagram is stored

    Returns:
        newVoronoi (vtkPolyData): Voronoi diagram
    """
    if path.isfile(filename):
        return read_polydata(filename)

    voronoi = vmtkscripts.vmtkDelaunayVoronoi()
    voronoi.Surface = surface
    voronoi.RemoveSubresolutionTetrahedra = 0
    voronoi.Execute()

    write_polydata(voronoi.VoronoiDiagram, filename)

    newVoronoi = voronoi.VoronoiDiagram
    return newVoronoi


def get_tolerance(centerline, N=50):
    """
    Finds tolerance based on
    average length between first N points
    along the input centerline.

    Args:
        centerline (vtkPolyData): Centerline data.
        N (int): Number of points

    Returns:
        tolerance (float): Tolerance value.
    """

    line = extract_single_line(centerline, 0)
    length = get_curvilinear_coordinate(line)
    tolerance = np.mean(length[1:N] - length[:N - 1]) / divergingRatioToSpacingTolerance

    return tolerance


def get_relevant_outlets(surface, dir_path):
    """
    Extract relevant outlets of the
    input surface model.

    Args:
        surface(vtkPolyData): Surface model.
        dir_path (str): Location of info-file.

    Returns:
        relevant_outlets (list): List of relevant outlet IDs.
    """
    # Check if info exists
    if not path.isfile(path.join(dir_path, "info.txt")):
        provide_relevant_outlets(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    relevant_outlets = []
    for key, value in list(parameters.items()):
        if key.startswith("relevant_outlet_"):
            relevant_outlets.append(value)

    if relevant_outlets == []:
        relevant_outlets = provide_relevant_outlets(surface, dir_path)

    return relevant_outlets


def smooth_voronoi_diagram(voronoi, centerlines, smoothingFactor,
                         no_smooth_cl=None):
    """
    Smooth voronoi diagram based on a given
    smoothingfactor. Each voronoi point
    that has a radius less then MISR*(1-smoothingFactor)
    at the closest centerline point is removed.

    Args:
        voronoi (vtkPolyData): Voronoi diagram to be smoothed.
        centerlines (vtkPolyData): Centerline data.
        smoothingFactor (float): Smoothing factor.
        no_smooth_cl (vktPolyData): Unsmoothed centerline.

    Returns: smoothedDiagram (vtkPolyData): Smoothed voronoi diagram.
    """
    numberOfPoints = voronoi.GetNumberOfPoints()

    threshold = get_array(radiusArrayName, centerlines) * (1 - smoothingFactor)
    locator = get_locator(centerlines)
    if no_smooth_cl is not None:
        no_locator = get_locator(no_smooth_cl)

    smoothedDiagram = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    radiusArrayNumpy = np.zeros(numberOfPoints)

    count = 0
    for i in range(numberOfPoints):
        point = voronoi.GetPoint(i)
        radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        id_ = locator.FindClosestPoint(point)
        cl_point = centerlines.GetPoint(id_)

        if distance(point, cl_point) > 3 * threshold[id_]:
            points.InsertNextPoint(point)
            cellArray.InsertNextCell(1)
            cellArray.InsertCellPoint(count)
            radiusArrayNumpy[count] = radius
            count += 1

        elif no_smooth_cl is not None:
            dist1 = distance(point, centerlines.GetPoint(id_))
            id_1 = no_locator.FindClosestPoint(point)
            dist2 = distance(point, no_smooth_cl.GetPoint(id_1))

            if dist2 < dist1:
                points.InsertNextPoint(point)
                cellArray.InsertNextCell(1)
                cellArray.InsertCellPoint(count)
                radiusArrayNumpy[count] = radius
                count += 1
            else:
                if radius >= threshold[id_]:
                    points.InsertNextPoint(point)
                    cellArray.InsertNextCell(1)
                    cellArray.InsertCellPoint(count)
                    radiusArrayNumpy[count] = radius
                    count += 1
        else:
            if radius >= threshold[id_]:
                points.InsertNextPoint(point)
                cellArray.InsertNextCell(1)
                cellArray.InsertCellPoint(count)
                radiusArrayNumpy[count] = radius
                count += 1

        radiusArray = get_vtk_array(radiusArrayName, 1, count)

    for i in range(count):
        radiusArray.SetTuple1(i, radiusArrayNumpy[i])

    smoothedDiagram.SetPoints(points)
    smoothedDiagram.SetVerts(cellArray)
    smoothedDiagram.GetPointData().AddArray(radiusArray)

    return smoothedDiagram


def get_curvilinear_coordinate(line):
    """
    Get curvilinear coordinates along
    an input centerline.

    Args:
        line (vtkPolyData): Input centerline

    Returns:
        curv_coor (ndarray): Array of abscissa points.
    """
    curv_coor = np.zeros(line.GetNumberOfPoints())
    for i in range(line.GetNumberOfPoints() - 1):
        pnt1 = np.asarray(line.GetPoints().GetPoint(i))
        pnt2 = np.asarray(line.GetPoints().GetPoint(i + 1))
        curv_coor[i + 1] = np.sqrt(np.sum((pnt1 - pnt2) ** 2)) + curv_coor[i]

    return curv_coor


def merge_data(inputs):
    """
    Appends one or more polygonal
    datates together into a single
    polygonal dataset.

    Args:
        inputs (list): List of vtkPolyData objects.

    Returns:
        merged_data (vtkPolyData): Single polygonal dataset.
    """
    appendFilter = vtk.vtkAppendPolyData()
    for input_ in inputs:
        appendFilter.AddInputData(input_)
    appendFilter.Update()
    merged_data = appendFilter.GetOutput()

    return merged_data


def get_array(arrayName, line, k=1):
    """
    Get data array from centerline object (Point data).

    Args:
        arrayName (str): Name of array.
        line (vtkPolyData): Centerline object.
        k (int): Dimension.

    Returns:
        array (ndarray): Array containing data points.
    """
    array = np.zeros((line.GetNumberOfPoints(), k))
    if k == 1:
        getData = line.GetPointData().GetArray(arrayName).GetTuple1
    elif k == 2:
        getData = line.GetPointData().GetArray(arrayName).GetTuple2
    elif k == 3:
        getData = line.GetPointData().GetArray(arrayName).GetTuple3

    for i in range(line.GetNumberOfPoints()):
        array[i, :] = getData(i)

    return array


def get_array_cell(arrayName, line, k=1):
    """
    Get data array from centerline object (Cell data).

    Args:
        arrayName (str): Name of array.
        line (vtkPolyData): Centerline object.
        k (int): Dimension.

    Returns:
        array (ndarray): Array containing data points.
    """
    array = np.zeros((line.GetNumberOfCells(), k))
    if k == 1:
        getData = line.GetCellData().GetArray(arrayName).GetTuple1
    elif k == 2:
        getData = line.GetCellData().GetArray(arrayName).GetTuple2
    elif k == 3:
        getData = line.GetCellData().GetArray(arrayName).GetTuple3
    elif k == 9:
        getData = line.GetCellData().GetArray(arrayName).GetTuple9

    for i in range(line.GetNumberOfCells()):
        array[i, :] = getData(i)

    return array


def create_new_surface(completeVoronoiDiagram, polyBallImageSize=[120, 120, 120]):
    """
    Envelops an input voronoi diagram
    into a new surface model at a
    given resolution determined by
    the polyBallImageSize.

    Args:
        completeVoronoiDiagram (vtkPolyData): Voronoi diagram
        polyBallImageSize (list): List of dimensional resolution of output model

    Returns:
        envelope (vtkPolyData): Enveloped surface model.
    """
    # TODO: Get bounds and compute polyballs based on that
    #       bounds = surface.GetBounds()
    # Get the x,y, and z range of the completeVoronoiDiagram
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    modeller.SetInputData(completeVoronoiDiagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(polyBallImageSize)
    modeller.Update()

    # Write the new surface
    marchingCube = vtk.vtkMarchingCubes()
    marchingCube.SetInputData(modeller.GetOutput())
    marchingCube.SetValue(0, 0.0)
    marchingCube.Update()
    envelope = marchingCube.GetOutput()

    return envelope


def get_aneurysm_dome(surface, dir_path, anu_num):
    """
    Extract aneurysm dome of the
    input surface model.

    Args:
        surface(vtkPolyData): Surface model.
        dir_path (str): Location of info-file.
        anu_num (int): Number of aneurysms.

    Returns:
        dome[anu_num] (list): List of aneurysm dome IDs.
    """
    # Check if info exists
    if not path.isfile(path.join(dir_path, "info.txt")):
        provide_aneurysm_points(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    dome = []
    for key, value in list(parameters.items()):
        if key.startswith("aneurysm_"):
            dome.append(value)

    # Recursive call
    if dome == []:
        dome = provide_aneurysm_points(surface, dir_path)

    return dome[anu_num]


def centerline_div(centerline1, centerline2, tol):
    """
    Find ID of diverging point;
    where two input centerlines diverge
    due to a bifurcation.

    Args:
        centerline1 (vtkPolyData): First centerline.
        centerline2 (vtkPolyData): Second centerline.
        tol (float): Tolerance.

    Returns:
        i (int): ID at diverging point.
    """
    # Find clipping points
    N = min(centerline1.GetNumberOfPoints(), centerline2.GetNumberOfPoints())
    get_point1 = centerline1.GetPoints().GetPoint
    get_point2 = centerline2.GetPoints().GetPoint

    for i in range(0, N):
        distance_between_points = distance(get_point1(i), get_point2(i))
        if distance_between_points > tol:
            break

    return i


def provide_relevant_outlets(surface, dir_path=None):
    """
    Get relevant outlets from user
    selected points on a input surface.

    Args:
        surface (vtkPolyData): Surface model.
        dir_path (str): Location of into.txt file

    Returns:
        points (list): List of relevant outlet IDs
    """

    # Fix surface
    cleaned_surface = surface_cleaner(surface)
    triangulated_surface = triangulate_surface(cleaned_surface)

    # Select seeds
    SeedSelector = vmtkPickPointSeedSelector()
    SeedSelector.SetSurface(triangulated_surface)
    SeedSelector.Execute()

    aneurysmSeedIds = SeedSelector.GetTargetSeedIds()
    get_point = surface.GetPoints().GetPoint
    points = [list(get_point(aneurysmSeedIds.GetId(i))) for i in range(aneurysmSeedIds.GetNumberOfIds())]
    info = {}

    if dir_path is not None:
        for i in range(len(points)):
            info["relevant_outlet_%d" % i] = points[i]
            write_parameters(info, dir_path)

    return points


def provide_aneurysm_points(surface, dir_path=None):
    """
    Get location of aneurysm(s) based on
    selected points on a input surface.

    Args:
        surface (vtkPolyData): Surface model.
        dir_path (str): Location of into.txt file

    Returns:
        points (list): List of aneurysm location IDs
    """
    # Fix surface
    cleaned_surface = surface_cleaner(surface)
    triangulated_surface = triangulate_surface(cleaned_surface)

    # Select seeds
    SeedSelector = vmtkPickPointSeedSelector()
    SeedSelector.SetSurface(triangulated_surface)
    SeedSelector.Execute()

    aneurysmSeedIds = SeedSelector.GetTargetSeedIds()
    get_point = surface.GetPoints().GetPoint
    points = [list(get_point(aneurysmSeedIds.GetId(i))) for i in range(aneurysmSeedIds.GetNumberOfIds())]

    if dir_path is not None:
        info = {"number_of_aneurysms": len(points)}

        for i in range(len(points)):
            info["aneurysm_%d" % i] = points[i]
            write_parameters(info, dir_path)

    return points


def get_data(centerline, centerline_bif, tol):
    """
    Locate bifurcating point and diverging points
    in a bifurcation.
    End points are set based on the MISR at
    the selected points.

    Args:
        centerline (vtkPolyData): Centerline from inlet to relevant outlets.
        centerline_bif (vtkPolyData): Centerline through bifurcation.
        tol (float): Tolerance parameter.

    Returns:
        data (dict): Contains info about diverging point locations.
    """
    # Sort centerline to start at inlet
    cl1 = extract_single_line(centerline, 0)
    cl2 = extract_single_line(centerline, 1)
    centerline = merge_data([cl1, cl2])

    # Declear variables before loop incase values are not found
    diverging_point_ID = -1
    diverging_point = [0.0, 0.0, 0.0]
    diverging_point_MISR = -1

    clipping_point_ID = -1
    clipping_point = [0.0, 0.0, 0.0]
    data = {"bif": {}, 0: {}, 1: {}}

    # List of points conected to ID
    points_ids_0 = vtk.vtkIdList()
    points_ids_1 = vtk.vtkIdList()

    # One is the branch to the left and the other is the one to the right
    centerline.GetCellPoints(0, points_ids_0)
    centerline.GetCellPoints(1, points_ids_1)
    # Find lower clipping point
    N = min(points_ids_0.GetNumberOfIds(), points_ids_1.GetNumberOfIds())
    for i in range(0, N):
        cell_point_0 = centerline.GetPoint(points_ids_0.GetId(i))
        cell_point_1 = centerline.GetPoint(points_ids_1.GetId(i))

        distance_between_points = distance(cell_point_0, cell_point_1) ** 2
        if distance_between_points > tol:
            tmpI = i
            point_ID_0 = points_ids_0.GetId(i)
            point_ID_1 = points_ids_1.GetId(i)
            center = centerline.GetPoint(point_ID_0)
            r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(point_ID_0)
            break

    end, r_end = move_past_sphere(centerline, center, r, point_ID_0, step=-1)
    data["bif"]["end_point"] = end
    data["bif"]["r_end"] = r_end
    data["bif"]["div_point"] = center
    data["bif"]["ID_div"] = point_ID_0
    data["bif"]["i_div"] = tmpI
    data["bif"]["r_div"] = r

    # Find the diverging points for the bifurcation
    # continue further downstream in each direction and stop when
    # a point is closer than tol, then move point MISR * X
    locator = get_locator(centerline_bif)

    for counter, point_ids in enumerate([points_ids_0, points_ids_1]):
        for i in range(tmpI, point_ids.GetNumberOfIds(), 1):
            tmp_point = centerline.GetPoint(point_ids.GetId(i))
            closest_point_ID = locator.FindClosestPoint(tmp_point)
            closest_point = centerline_bif.GetPoint(closest_point_ID)
            distance_between_points = distance(tmp_point, closest_point) ** 2
            if distance_between_points < tol:
                point_ID = point_ids.GetId(i)
                center = centerline.GetPoint(point_ID)
                r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(point_ID)
                break

        end, r_end = move_past_sphere(centerline, center, r, point_ID, step=1, stop=point_ID * 100, X=1)
        data[counter]["end_point"] = end
        data[counter]["r_end"] = r_end
        data[counter]["r_div"] = r
        data[counter]["ID_end"] = locator.FindClosestPoint(data[counter]["end_point"])
        data[counter]["ID_div"] = locator.FindClosestPoint(center)
        data[counter]["div_point"] = center

    return data


def write_points(points, filename):
    """
    Writes input points to file.

    Args:
        points (vtkPolyData): Point data.
        filename (str): Save location.
    """
    pointSet = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()

    for i in range(points.GetNumberOfPoints()):
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i)

    pointSet.SetPoints(points)
    pointSet.SetVerts(cellArray)

    write_polydata(pointSet, filename)


def surface_cleaner(surface):
    """
    Clean surface by merging
    duplicate points, and∕or
    removing unused points
    and∕or removing degenerate cells.

    Args:
        surface (vtkPolyData): Surface model.

    Returns:
        cleanSurface (vtkPolyData): Cleaned surface model.
    """
    # Clean surface
    surfaceCleaner = vtk.vtkCleanPolyData()
    surfaceCleaner.SetInputData(surface)
    surfaceCleaner.Update()
    cleanSurface = surfaceCleaner.GetOutput()

    return cleanSurface


def get_centers(surface, dir_path, flowext=False):
    # Check if info exists
    if flowext or not path.isfile(path.join(dir_path, "info.txt")):
        compute_centers(surface, dir_path)

    # Open info
    parameters = get_parameters(dir_path)
    outlets = []
    inlet = []
    for key, value in list(parameters.items()):
        if key == "inlet":
            inlet = value
        elif "outlet" in key and "area" not in key and "relevant" not in key:
            outlets += value

    num_outlets = len(outlets) // 3
    if num_outlets != 0:
        outlets = []
        for i in range(num_outlets):
            outlets += parameters["outlet%d" % i]

    if inlet == [] and outlets == []:
        inlet, outlets = compute_centers(surface, dir_path)

    return inlet, outlets


def triangulate_surface(surface):
    """Triangulate surface or polygon(?)"""
    surfaceTriangulator = vtk.vtkTriangleFilter()
    surfaceTriangulator.SetInputData(surface)
    surfaceTriangulator.PassLinesOff()
    surfaceTriangulator.PassVertsOff()
    surfaceTriangulator.Update()

    return surfaceTriangulator.GetOutput()


def geometry_filter(unstructured_grid):
    # Convert unstructured grid to polydata
    filter = vtk.vtkGeometryFilter()
    filter.SetInputData(unstructured_grid)
    filter.Update()
    polydata = filter.GetOutput()

    return polydata


def threshold(surface, name, lower=0, upper=1, type="between", source=1):
    # source = 1 uses cell data as input
    # source = 0 uses point data as input

    # Apply threshold
    threshold = vtk.vtkThreshold()
    threshold.SetInputData(surface)
    if type == "between":
        threshold.ThresholdBetween(lower, upper)
    elif type == "lower":
        threshold.ThresholdByLower(lower)
    elif type == "upper":
        threshold.ThresholdByUpper(upper)
    else:
        print((("%s is not a threshold type. Pleace chose from: upper, lower" + \
              ", or between") % type))
        sys.exit(0)

    threshold.SetInputArrayToProcess(0, 0, 0, source, name)
    threshold.Update()
    surface = threshold.GetOutput()

    # Convert to polydata
    surface = geometry_filter(surface)

    return surface


def compute_area(surface):
    "Compute area of polydata"
    mass = vtk.vtkMassProperties()
    mass.SetInputData(surface)

    return mass.GetSurfaceArea()


def uncapp_surface(surface):
    # Add-hoc method for removing capps on surfaces
    # This could proboly be highly improved, but is sufficient for now.

    # Get cell normals
    normal_generator = vtk.vtkPolyDataNormals()
    normal_generator.SetInputData(surface)
    normal_generator.ComputePointNormalsOff()
    normal_generator.ComputeCellNormalsOn()
    normal_generator.Update()
    cell_normals = normal_generator.GetOutput()

    # Compute gradients of the normals
    gradients_generator = vtk.vtkGradientFilter()
    gradients_generator.SetInputData(cell_normals)
    gradients_generator.SetInputArrayToProcess(0, 0, 0, 1, "Normals")
    gradients_generator.Update()
    gradients = gradients_generator.GetOutput()

    # Compute the magnitude of the gradient
    gradients_array = get_array_cell("Gradients", gradients, k=9)
    gradients_magnitude = np.sqrt(np.sum(gradients_array ** 2, axis=1))

    # Mark all cells with a gradient magnitude less then 0.1
    end_capp_array = gradients_magnitude < 0.08
    end_capp_vtk = get_vtk_array("Gradients_mag", 1, end_capp_array.shape[0])
    for i, p in enumerate(end_capp_array):
        end_capp_vtk.SetTuple(i, [p])
    gradients.GetCellData().AddArray(end_capp_vtk)

    # Extract capps
    end_capps = threshold(gradients, "Gradients_mag", lower=0.5, upper=1.5,
                            type="between", source=1)

    # Get connectivity
    end_capps_connectivity = get_connectivity(end_capps)
    region_array = get_array("RegionId", end_capps_connectivity)

    # Compute area for each region
    area = []
    circleness = []
    regions = []
    centers_edge = []
    limit = 0.1
    for i in range(int(region_array.max()) + 1):
        regions.append(threshold(end_capps_connectivity, "RegionId", lower=(i - limit),
                            upper=(i + limit), type="between", source=0))
        circ, center = compute_circleness(regions[-1])
        circleness.append(circ)
        centers_edge.append(center)
        area.append(compute_area(regions[-1]))

    # Only keep outlets with circleness < 3 and area > 0.3 mm^2
    circleness_IDs = np.where(np.array(circleness) < 3)
    region_IDs = np.where(np.array(area) > 0.3)
    regions = [regions[i] for i in region_IDs[0] if i in circleness_IDs[0]]
    centers_edge = [centers_edge[i] for i in region_IDs[0] if i in circleness_IDs[0]]

    # Mark the outlets on the original surface
    mark_outlets = create_vtk_array(np.zeros(surface.GetNumberOfCells()), "outlets", k=1)
    locator = get_locator_cell(surface)
    tmp_center = [0, 0, 0]
    for region in regions:
        centers_filter = vtk.vtkCellCenters()
        centers_filter.SetInputData(region)
        centers_filter.VertexCellsOn()
        centers_filter.Update()
        centers = centers_filter.GetOutput()

        for i in range(centers.GetNumberOfPoints()):
            centers.GetPoint(i, tmp_center)
            p = [0, 0, 0]
            cellId = vtk.mutable(0)
            subId = vtk.mutable(0)
            dist = vtk.mutable(0)
            locator.FindClosestPoint(tmp_center, p, cellId, subId, dist)
            mark_outlets.SetTuple(cellId, [1])

    surface.GetCellData().AddArray(mark_outlets)

    # Remove the outlets from the original surface
    uncapped_surface = threshold(surface, "outlets", lower=0, upper=0.5,
                                type="between", source=1)

    # Check if some cells where not marked
    remove = True
    while remove:
        locator = get_locator_cell(uncapped_surface)
        mark_outlets = create_vtk_array(np.zeros(uncapped_surface.GetNumberOfCells()), "outlets", k=1)
        remove = False
        for center in centers_edge:
            locator.FindClosestPoint(center, p, cellId, subId, dist)
            if dist < 0.01:
                remove = True
                mark_outlets.SetTuple(cellId, [1])

        uncapped_surface.GetCellData().AddArray(mark_outlets)

        if remove:
            uncapped_surface = threshold(uncapped_surface, "outlets", lower=0,
                                            upper=0.5, type="between", source=1)

    return uncapped_surface


def capp_surface(surface):
    surfaceCapper = vtkvmtk.vtkvmtkCapPolyData()
    surfaceCapper.SetInputData(surface)
    surfaceCapper.SetDisplacement(0.0)
    surfaceCapper.SetInPlaneDisplacement(0.0)
    surfaceCapper.Update()

    return surfaceCapper.GetOutput()


def is_surface_capped(surface):
    return compute_centers(surface, test_capped=True)


def get_connectivity(surface, mode="All", closestPoint=None):
    """Compute connectivity of the cells"""
    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(surface)

    # Mark each region with "RegionId"
    if mode == "All":
        connectivity.SetExtractionModeToAllRegions()
    elif mode == "Largest":
        connectivity.SetExtractionModeToLargestRegion()
    elif mode == "Closest":
        if closestPoint is None:
            print("ERROR: point not set for extracting closest region")
            sys.exit(0)
        connectivity.SetExtractionModeToClosestPointRegion()
        connectivity.SetClosestPoint(closestPoint)
    connectivity.ColorRegionsOn()
    connectivity.Update()
    output = connectivity.GetOutput()
    return output


def compute_circleness(surface):
    edges = getFeatureEdges(surface)

    # Get points
    points = []
    for i in range(edges.GetNumberOfPoints()):
        points.append(edges.GetPoint(i))

    # Compute center
    points = np.array(points)
    center = np.mean(np.array(points), axis=0)

    # Compute ratio between max inscribed sphere, and min inscribed "area"
    point_radius = np.sqrt(np.sum((points - center) ** 2, axis=1))
    argsort = np.argsort(point_radius)
    if point_radius[argsort[1]] / point_radius[argsort[0]] > 15:
        radius_min = point_radius[argsort[1]]
    else:
        radius_min = point_radius.min()

    min_area = math.pi * radius_min ** 2
    max_area = math.pi * point_radius.max() ** 2

    return max_area / min_area, center


def get_feature_edges(polyData):
    """Extracts the edges of the cells that are open"""
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOff()
    featureEdges.SetInputData(polyData)
    featureEdges.Update()

    return featureEdges.GetOutput()


def compute_centers(polyData, case_path=None, test_capped=False):
    # Get cells which are open
    cells = get_feature_edges(polyData)

    # Check is the model is closed
    if test_capped:
        if cells.GetNumberOfCells() == 0:
            return True
        else:
            return False

    if cells.GetNumberOfCells() == 0:
        print("WARNING: The model is capped, so it is uncapped, but the method is experimental.")
        uncapped_surface = uncapp_surface(polyData)
        compute_centers(uncapped_surface, case_path, test_capped)

    # Compute connectivity of the cells
    outputs = get_connectivity(cells)

    # Get connectivity array
    region_array = get_array("RegionId", outputs)

    if test_capped:
        return region_array.max() >= 1

    # Get points
    points = []
    get_point = outputs.GetPoints().GetPoint
    for i in range(region_array.shape[0]):
        points.append(get_point(i))
    points = np.asarray(points)

    # Get area and center
    area = []
    center = []
    for i in range(int(region_array.max()) + 1):
        # Extract points for this opening
        index = (region_array == i).nonzero()[0]
        n = index.shape[0]
        tmp_points = np.zeros((n, 3))
        if n <= 5:
            continue

        # Sort to be ordred clockwise
        for j in range(3):
            x = points[index][:, j]
            x_end = x[2::2].tolist()
            x_end.reverse()
            x1 = [x[0]] + x[1::2].tolist() + x_end
            tmp_points[:, j] = x1

        # Insert points into VTK object
        points_vtk = vtk.vtkPoints()
        [points_vtk.InsertNextPoint(tmp_points[k]) for k in range(n)]

        # Make polygon
        polygon = vtk.vtkPolygon()
        polygon.GetPoints().DeepCopy(points_vtk)
        polygon.GetPointIds().SetNumberOfIds(n)
        for j in range(n):
            polygon.GetPointIds().SetId(j, j)

        area.append(polygon.ComputeArea())
        # FIXME: Remve bary center as well
        center.append(np.array(compute_bary_center(tmp_points)))

    # Store the center and area
    inlet_ind = area.index(max(area))
    if case_path is not None:
        info = {"inlet": center[inlet_ind].tolist(), "inlet_area": area[inlet_ind]}
        p = 0
        for i in range(len(area)):
            if i == inlet_ind:
                p = -1
                continue
            info["outlet%d" % (i + p)] = center[i].tolist()
            info["outlet%s_area" % (i + p)] = area[i]

        write_parameters(info, case_path)

    inlet_center = center[inlet_ind].tolist()
    center.pop(inlet_ind)

    center_ = [item for sublist in center for item in sublist]

    return inlet_center, center_


def compute_bary_center(points):
    # Get i+1
    shifted = np.zeros(points.shape)
    shifted[1:, :] = points[:-1, :]
    shifted[0, :] = points[-1, :]

    # Compute weights
    weight = np.sqrt(np.sum((points - shifted) ** 2, axis=1))
    weight_sum = np.sum(weight)

    # Compute center
    center_x = np.sum((points[:, 0] + shifted[:, 0]) / 2 * weight) / weight_sum
    center_y = np.sum((points[:, 1] + shifted[:, 1]) / 2 * weight) / weight_sum
    center_z = np.sum((points[:, 2] + shifted[:, 2]) / 2 * weight) / weight_sum

    return [center_x, center_y, center_z]


def get_vtk_array(name, comp, num):
    array = vtk.vtkDoubleArray()
    array.SetNumberOfComponents(comp)
    array.SetNumberOfTuples(num)
    for i in range(comp):
        array.FillComponent(i, 0.0)
    array.SetName(name)

    return array


def get_locator_cell(surface):
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()

    return locator


def get_locator(centerline):
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(centerline)
    locator.BuildLocator()

    return locator


def distance(point1, point2):
    return np.sqrt(np.sum((np.asarray(point1) - np.asarray(point2)) ** 2))


def remove_distant_points(voronoi, centerline, limit=None):
    N = voronoi.GetNumberOfPoints()
    newVoronoi = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    radius = np.zeros(N)

    locator = get_locator(centerline)
    get_data = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1

    if limit is None:
        limit = []
        m = N // 500
        steps = list(range(N))[::m]
        for i in steps:
            if 1 < get_data(i) < 10:
                limit.append(get_data(i))
        limit_ = np.median(limit) * 3

    count = 0
    for i in range(N):
        point = voronoi.GetPoint(i)
        ID = locator.FindClosestPoint(point)
        cl_point = centerline.GetPoint(ID)
        dist = distance(point, cl_point)
        if limit(get_data(i), point, dist):
            count += 1
            continue

        points.InsertNextPoint(point)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i - count)
        value = get_data(i)
        radius[i - count] = value

    print("Removed %s points from the voronoi diagram" % count)

    radiusArray = get_vtk_array(radiusArrayName, 1, N - count)
    for i in range(N - count):
        radiusArray.SetTuple(i, [float(radius[i])])

    newVoronoi.SetPoints(points)
    newVoronoi.SetVerts(cellArray)
    newVoronoi.GetPointData().AddArray(radiusArray)

    return newVoronoi


def vmtk_compute_centerline_sections(surface, centerline):
    centerlineSections = vtkvmtk.vtkvmtkPolyDataCenterlineSections()
    centerlineSections.SetInputData(surface)
    centerlineSections.SetCenterlines(centerline)
    centerlineSections.SetCenterlineSectionAreaArrayName('CenterlineSectionArea')
    centerlineSections.SetCenterlineSectionMinSizeArrayName('CenterlineSectionMinSize')
    centerlineSections.SetCenterlineSectionMaxSizeArrayName('CenterlineSectionMaxSize')
    centerlineSections.SetCenterlineSectionShapeArrayName('CenterlineSectionShape')
    centerlineSections.SetCenterlineSectionClosedArrayName('CenterlineSectionClosed')
    centerlineSections.Update()

    CenterlineSections = centerlineSections.GetOutput()
    line = centerlineSections.GetCenterlines()

    return line, CenterlineSections


def vmtk_compute_centerlines(inlet, outlet, filepath, surface, resampling=1,
                        smooth=False, num_iter=100, smooth_factor=0.1,
                        endPoint=1, method="pointlist"):
    if path.isfile(filepath):
        return read_polydata(filepath)

    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.SeedSelectorName = 'pointlist'
    centerlines.AppendEndPoints = endPoint
    centerlines.Resampling = 1
    centerlines.ResamplingStepLength = resampling
    centerlines.SourcePoints = inlet
    centerlines.TargetPoints = outlet
    centerlines.Execute()
    centerlines = centerlines.Centerlines

    if smooth:
        centerlineSmoothing = vmtkscripts.vmtkCenterlineSmoothing()
        centerlineSmoothing.SetInputData(self.Centerlines)
        centerlineSmoothing.SetNumberOfSmoothingIterations(num_iter)
        centerlineSmoothing.SetSmoothingFactor(smooth_factor)
        centerlineSmoothing.Update()

        centerlines = centerlinesSmooth.GetOutput()

    # Save the computed centerline.
    if filepath is not None:
        write_polydata(centerlines, filepath)

    return centerlines


def create_vtk_array(values, name, k=1):
    vtkArray = get_vtk_array(name, k, values.shape[0])

    if k == 1:
        for i in range(values.shape[0]):
            vtkArray.SetTuple1(i, values[i])
    elif k == 2:
        for i in range(values.shape[0]):
            vtkArray.SetTuple2(i, values[i, 0], values[i, 1])
    elif k == 3:
        for i in range(values.shape[0]):
            vtkArray.SetTuple3(i, values[i, 0], values[i, 1], values[i, 2])
    elif k == 9:
        for i in range(values.shape[0]):
            vtkArray.SetTuple9(i, values[i, 0], values[i, 1], values[i, 2],
                                  values[i, 3], values[i, 4], values[i, 5],
                                  values[i, 6], values[i, 7], values[i, 8])

    return vtkArray


def gram_schmidt(V):
    V = 1.0 * V
    U = np.copy(V)

    def proj(u, v):
        return u * np.dot(v, u) / np.dot(u, u)

    for i in range(1, V.shape[1]):
        for j in range(i):
            U[:, i] -= proj(U[:, j], V[:, i])

    # normalize column
    den = (U ** 2).sum(axis=0) ** 0.5
    E = U / den

    return E


def get_parameters(folder):
    # If info.txt file, return an empty dict
    if not path.isfile(path.join(folder, "info.txt")):
        return {}

    # Get text
    f = open(path.join(folder, "info.txt"), "r")
    text = f.read()
    f.close()
    text = text.split("\n")

    # Get values
    data = {}
    for par in text:
        if par != "":
            split = par.split(": ")
            if len(split) == 2:
                key, value = split
            else:
                key = split[0]
                value = ": ".join(split[1:])
            try:
                data[key] = eval(value)
            except:
                data[key] = value
    return data


def write_parameters(data, folder):
    # Get old parameters
    parameters = get_parameters(folder)

    # Add new parameters (can overwrite old as well)
    for key, value in list(data.items()):
        parameters[key] = value

    # Get new text
    text = ["%s: %s" % (k, v) for k, v in list(parameters.items())]
    text = "\n".join(text)

    # Write text
    f = open(path.join(folder, "info.txt"), "w")
    f.write(text)
    f.close()


def data_to_vtkPolyData(data, header, TNB=None, PT=None):
    line = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    cellArray.InsertNextCell(data.shape[0])
    linePoints = vtk.vtkPoints()

    info_array = []
    for i in range(3, data.shape[1]):
        radiusArray = get_vtk_array(header[i], 1, data.shape[0])
        info_array.append(radiusArray)

    if TNB is not None:
        for i in range(3):
            radiusArray = get_vtk_array(header[i + data.shape[1]], 3, data.shape[0])
            info_array.append(radiusArray)

    if PT is not None:
        start = data.shape[1] if TNB is None else data.shape[1] + 3
        for i in range(2):
            radiusArray = get_vtk_array(header[i + start], 3, PT[0].shape[0])
            info_array.append(radiusArray)

    for i in range(data.shape[0]):
        cellArray.InsertCellPoint(i)
        linePoints.InsertNextPoint(data[i, :3])
        for j in range(3, data.shape[1]):
            info_array[j - 3].SetTuple1(i, data[i, j])

    if TNB is not None:
        for i in range(data.shape[0]):
            for j in range(data.shape[1] - 3, data.shape[1], 1):
                tnb_ = TNB[j - data.shape[1]][i, :]
                info_array[j].SetTuple3(i, tnb_[0], tnb_[1], tnb_[2])

    if PT is not None:
        start = data.shape[1] - 3 if TNB is None else data.shape[1]
        for i in range(PT[-1].shape[0]):
            for j in range(start, start + 2, 1):
                pt_ = PT[j - start][i, :]
                info_array[j].SetTuple3(i, pt_[0], pt_[1], pt_[2])

    line.SetPoints(linePoints)
    line.SetLines(cellArray)
    for i in range(len(header) - 3):
        line.GetPointData().AddArray(info_array[i])

    return line


def get_number_of_arrays(line):
    count = 0
    names = []
    name = 0
    while name is not None:
        name = line.GetPointData().GetArrayName(count)
        if name is not None:
            names.append(name)
            count += 1

    return count, names


def extract_single_line(centerlines, id, startID=0, endID=None):
    cell = vtk.vtkGenericCell()
    centerlines.GetCell(id, cell)
    N = cell.GetNumberOfPoints() if endID is None else endID + 1

    line = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    cellArray.InsertNextCell(N - startID)
    linePoints = vtk.vtkPoints()

    arrays = []
    N_, names = get_number_of_arrays(centerlines)
    for i in range(N_):
        tmp = centerlines.GetPointData().GetArray(names[i])
        tmp_comp = tmp.GetNumberOfComponents()
        radiusArray = get_vtk_array(names[i], tmp_comp, N - startID)
        arrays.append(radiusArray)

    getArray = []
    for i in range(N_):
        getArray.append(centerlines.GetPointData().GetArray(names[i]))

    count = 0
    for i in range(startID, N):
        cellArray.InsertCellPoint(count)
        linePoints.InsertNextPoint(cell.GetPoints().GetPoint(i))

        for j in range(N_):
            num = getArray[j].GetNumberOfComponents()
            if num == 1:
                tmp = getArray[j].GetTuple1(i)
                arrays[j].SetTuple1(count, tmp)
            elif num == 2:
                tmp = getArray[j].GetTuple2(i)
                arrays[j].SetTuple2(count, tmp[0], tmp[1])
            elif num == 3:
                tmp = getArray[j].GetTuple3(i)
                arrays[j].SetTuple3(count, tmp[0], tmp[1], tmp[2])
            elif num == 9:
                tmp = getArray[j].GetTuple9(i)
                arrays[j].SetTuple9(count, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4],
                                       tmp[5], tmp[6], tmp[7], tmp[8])
        count += 1

    line.SetPoints(linePoints)
    line.SetLines(cellArray)
    for j in range(N_):
        line.GetPointData().AddArray(arrays[j])

    return line


def move_past_sphere(centerline, center, r, start, step=-1, stop=0, X=0.8):
    """Moves a point along the centerline until it as outside MIS"""
    # Create the minimal inscribed sphere
    MISphere = vtk.vtkSphere()
    MISphere.SetCenter(center)
    MISphere.SetRadius(r * X)
    tempPoint = [0.0, 0.0, 0.0]

    # Go the length of one MISR backwards
    for i in range(start, stop, step):
        value = MISphere.EvaluateFunction(centerline.GetPoint(i))
        if value >= 0.0:
            tempPoint = centerline.GetPoint(i)
            break

    r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

    return tempPoint, r


def vmtk_surface_smoother(surface, method, iterations=800):
    smoother = vmtkscripts.vmtkSurfaceSmoothing()
    smoother.Surface = surface
    smoother.NumberOfIterations = iterations
    smoother.Method = method
    smoother.Execute()
    surface = smoother.Surface

    return surface


def make_centerline(ifile, ofile, length=0.1, it=100, factor=0.1, in_out=None,
                   smooth=False, resampling=False, newpoints=None,
                   recompute=False, store_points=False, endpoints=0):
    """
    A general centerline command. If a centerline file with the same file
    name alread exists, then the file is just read. To overwrite this set
    recompute to True. If recomputed is to True then it uses the exsisting points
    from the old centerline file and no interaction with the interface is
    needed. Further on one can choose witch points you want to include by
    giving in_out. The first element in the list is the source point, and if -1
    is given, is chooses the old inlet, else it chooses outletX, where X is the
    outlet ID.

    Args:
        ifile (vtkPolyData): Input surface model.
        ofile (str): Name of output file.
        length (float): Number in (0,1], the resampling step.
        it(int):  Number of iterations of smoothing
        factor (float): The smoothening factor
        smooth (bool): Used to turn on/off smoothening
        resampling (bool): Used to turn in/off resampling
        in_out (list): Outlet/inlet points that should be used when recompute
        newpoints (list): Outlet/inlet points that should be used when recompute
    """

    # Check if ifile exsists
    if not path.exists(ifile):
        print("The input file: %s does not exsists!" % ifile)
        sys.exit(0)

    # Check if it already exsists or if it is to be recomputed
    if not path.exists(ofile) or recompute:
        # If recomputed use the old source and target points
        basedir = path.sep.join(path.join(ifile.split(path.sep)[:-2]))
        parameters = get_parameters(basedir)

        vmtkcenterlines = vmtkscripts.vmtkCenterlines()
        vmtkcenterlines.Surface = read_polydata(ifile)

        if newpoints is not None:
            inlet = newpoints[-1]
            points_ = newpoints[:-1]

        elif "inlet" in parameters:
            if in_out is None:
                inlet = parameters["inlet"]
                out = [k for k in list(parameters.keys()) if "outlet" in k and len(k) < 12]
                out.sort()
                points_ = [parameters[p] for p in out]
            else:
                inlet = parameters["inlet"] if in_out[0] == -1 else parameters["outlet%s" % in_out[0]]
                points_ = [parameters["outlet%s" % i] for i in in_out[1:]]

        if newpoints is not None or "inlet" in parameters:
            # Extract outlet points
            outlets = []
            for p in points_:
                for p_ in p:
                    outlets.append(p_)

            vmtkcenterlines.SeedSelectorName = "pointlist"
            vmtkcenterlines.SourcePoints = inlet
            vmtkcenterlines.TargetPoints = outlets
        else:
            vmtkcenterlines.SeedSelectorName = "pickpoint"

        # Add resampling
        if resampling:
            vmtkcenterlines.Resampling = 1
            vmtkcenterlines.ResamplingStepLength = length

        # Exectue command and save centerline
        vmtkcenterlines.AppendEndPoints = endpoints
        vmtkcenterlines.Execute()

        centerline = vmtkcenterlines.Centerlines
        write_polydata(centerline, ofile)

        # Add smoothing
        if smooth:
            centerlineSmoothing = vmtkscripts.vmtkCenterlineSmoothing()
            centerlineSmoothing.Centerlines = centerline
            centerlineSmoothing.NumberOfSmoothingIterations = it
            centerlineSmoothing.SmoothingFactor = factor
            centerlineSmoothing.Execute()
            centerline = centerlinesSmooth.Centerlines

        # If the points are not already stored, do it now
        if store_points:
            for i in range(centerline.GetNumberOfLines()):
                tmp_line = extract_single_line(centerline, i)
                tmp_N = tmp_line.GetNumberOfPoints()
                parameters["outlet%s" % i] = tmp_line.GetPoint(tmp_N - 1)
            parameters["inlet"] = tmp_line.GetPoint(0)
            write_parameters(parameters, basedir)
    else:
        centerline = read_polydata(ofile)

    return centerline


# Extract carotid siphon
def extract_carotid_siphon(folder):
    centerline_path = path.join(folder, "surface", "model_usr_centerline.vtp")
    centerline_bif_path = path.join(folder, "surface", "model_usr_centerline_bif.vtp")
    centerline_bif = read_polydata(centerline_bif_path)
    centerline = read_polydata(centerline_path)

    centerlineSpacing = np.sqrt(distance(centerline.GetPoint(2), centerline.GetPoint(3)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance
    data = get_data(centerline, centerline_bif, centerlineSpacing)
    line = extract_single_line(centerline, 0, startID=0, endID=data["bif"]["ID_div"])
    write_polydata(line, path.join(folder, "surface", "carotid_siphon.vtp"))

    return line


# Resampling of centerline
def vmtk_centerline_resampling(line, length, filename=None):
    resampler = vmtkscripts.vmtkCenterlineResampling()
    resampler.Centerlines = line
    resampler.Length = length
    resampler.Execute()

    line = resampler.Centerlines

    return line


# Check bool value for argparse
def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


# Compute geometric features of centerline
def vmtk_centerline_geometry(line, smooth, outputsmoothed=False, factor=1.0, iterations=100):
    geometry = vmtkscripts.vmtkCenterlineGeometry()
    geometry.Centerlines = line
    if smooth:
        geometry.LineSmoothing = 1
        geometry.OutputSmoothedLines = outputsmoothed
        geometry.SmoothingFactor = factor
        geometry.NumberOfSmoothingIterations = iterations
    geometry.CurvatureArrayName = "Curvature"
    geometry.TorsionArrayName = "Torsion"
    geometry.Execute()

    line = geometry.Centerlines
    return line


def vmtk_centerline_attributes(line):
    attributes = vmtkscripts.vmtkCenterlineAttributes()
    attributes.Centerlines = line
    attributes.NormalsArrayName = parallelTransportNormalsArrayName
    attributes.AbscissaArrayName = abscissasArrayName
    attributes.Execute()

    line = attributes.Centerlines
    return line


def discrete_geometry(line, neigh=10):
    len_line = line.GetNumberOfPoints()
    N = line.GetNumberOfPoints()

    # Compute cumulative chord length
    t = np.zeros(N)
    p = []
    for i in range(N):
        p.append(np.array(list(line.GetPoint(i))))
        p[i] = np.array(p[i])

    norms = [la.norm(p[j] - p[j - 1]) for j in range(1, N)]
    s = sum(norms)
    for i in range(1, N):
        s1 = sum(norms[:i + 1])
        t[i] = s1 / s

    # Radius of sliding neighbourhood
    m = neigh

    dxdt = np.zeros(N)
    dydt = np.zeros(N)
    dzdt = np.zeros(N)

    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)

    for i in range(N):
        x[i] = p[i][0]
        y[i] = p[i][1]
        z[i] = p[i][2]

    for i in range(0, m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(0, 2 * m + 1)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(0, 2 * m + 1)]) / t_sum

    for i in range(m, N - m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(i - m, i + m + 1)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(i - m, i + m + 1)]) / t_sum

    for i in range(N - m, N):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(N - 2 * m, N)])
        dxdt[i] = sum([(t[j] - t[i]) * (x[j] - x[i]) for j in range(N - 2 * m - 1, N)]) / t_sum
        dydt[i] = sum([(t[j] - t[i]) * (y[j] - y[i]) for j in range(N - 2 * m - 1, N)]) / t_sum
        dzdt[i] = sum([(t[j] - t[i]) * (z[j] - z[i]) for j in range(N - 2 * m - 1, N)]) / t_sum

    dgammadt = []
    dgammadt_norm = np.zeros(N)
    for i in range(N):
        dgammadt.append(np.array([dxdt[i], dydt[i], dzdt[i]]))
        dgammadt_norm[i] = la.norm(dgammadt[i])

    tg = []
    for i in range(N):
        tg.append(dgammadt[i] / dgammadt_norm[i])

    t1 = np.zeros(N)
    t2 = np.zeros(N)
    t3 = np.zeros(N)

    for i in range(N):
        t1[i] = tg[i][0]
        t2[i] = tg[i][1]
        t3[i] = tg[i][2]

    dt1dt = np.zeros(N)
    dt2dt = np.zeros(N)
    dt3dt = np.zeros(N)

    for i in range(0, m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(0, 2 * m + 1)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(0, 2 * m + 1)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(0, 2 * m + 1)]) / t_sum

    for i in range(m, N - m):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(i - m, i + m + 1)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(i - m, i + m + 1)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(i - m, i + m + 1)]) / t_sum

    for i in range(N - m, N):
        t_sum = sum([(t[j] - t[i]) ** 2 for j in range(N - 2 * m, N)])
        dt1dt[i] = sum([(t[j] - t[i]) * (t1[j] - t1[i]) for j in range(N - 2 * m - 1, N)]) / t_sum
        dt2dt[i] = sum([(t[j] - t[i]) * (t2[j] - t2[i]) for j in range(N - 2 * m - 1, N)]) / t_sum
        dt3dt[i] = sum([(t[j] - t[i]) * (t3[j] - t3[i]) for j in range(N - 2 * m - 1, N)]) / t_sum

    dtgdt = []
    dtgdt_norm = np.zeros(N)
    for i in range(N):
        dtgdt.append(np.array([dt1dt[i], dt2dt[i], dt3dt[i]]))
        dtgdt_norm[i] = la.norm(dtgdt[i])

    curv = np.zeros(N)
    for i in range(N):
        curv[i] = dtgdt_norm[i] / dgammadt_norm[i]

    curv = resample(curv, len_line)

    return line, curv


def get_k1k2_basis(curvature, line):
    """
    Create a k1-k2 basis used to determine
    the location of each bend of the carotid
    siphon.

    Args:
        curvature (floats): Curvature array.
        line (vtk): Centerline points.

    Returns:
        line (vtk): Centerline points including k1-k2 basis.
    """

    # The ParallelTransportNormals and the FrenetTangent is not orthonormal
    # (but close) from vmtk. Using the GramSchmidt proses gives E2 and fixes
    # the non-orthogonality
    E1 = get_array("ParallelTransportNormals", line, k=3)
    T = get_array("FrenetTangent", line, k=3)
    E2 = np.zeros((E1.shape[0], 3))

    for i in range(E1.shape[0]):
        V = np.eye(3)
        V[:, 0] = T[i, :]
        V[:, 1] = E1[i, :]
        V = gram_schmidt(V)

        E1[i, :] = V[:, 1]
        E2[i, :] = V[:, 2]

    # Compute k_1, k_2 furfilling T' = curv(s)*N(s) = k_1(s)*E_1(s) + k_2(s)*E_2(s).
    # This is simply a change of basis for the curvature vector N. The room
    # of k_1 and k_2 can be used to express both the curvature and the
    # torsion.
    N = get_array("FrenetNormal", line, k=3)

    k2 = (curvature.T * (E1[:, 1] * N[:, 0] - N[:, 1] * E1[:, 0]) / \
                        (E2[:, 1] * E1[:, 0] - E2[:, 0] * E1[:, 1]))[0]
    k1 = (-(curvature.T * N[:, 0] + k2 * E2[:, 0]) / E1[:, 0])[0]

    for k in [(k1, "k1"), (k2, "k2")]:
        k_array = create_vtk_array(k[0], k[1])
        line.GetPointData().AddArray(k_array)

    return line


def spline_centerline(line, get_curv=False, isline=False, nknots=50, get_stats=True):
    """
    Given the knots and coefficients of a B-spline representation,
    evaluate the value of the smoothing polynomial and its derivatives.
    This is a wrapper around the FORTRAN routines splev and splder of FITPACK.

    Args:
        line (vtkPolyData): Centerline points.
        get_curv (bool): Computes curvature profile if True.
        isline (bool): Determines if centerline object is a line or points.
        nknots (int): Number of knots.
        get_stats (bool): Determines if curve attribuites are computed or not.

    Returns:
        line (vtkPolyData): Splined centerline data.
    Returns:
        curv (ndarray): Curvature profile.
    """

    if not isline:
        # Create edges between new_centerline points
        pts = vtk.vtkPoints()
        for i in range((line.GetNumberOfCells())):
            pts.InsertNextPoint(line.GetPoint(i))

        lines = vtk.vtkCellArray()
        for i in range((line.GetNumberOfCells()) - 1):
            newline = vtk.vtkLine()
            newline.GetPointIds().SetId(0, i)
            newline.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(newline)

        line = vtk.vtkPolyData()
        line.SetPoints(pts)
        line.SetLines(lines)

    # Collect data from centerline
    data = np.zeros((line.GetNumberOfPoints(), 3))
    MISR = get_array(radiusArrayName, line)

    curv_coor = get_curvilinear_coordinate(line)
    for i in range(data.shape[0]):
        data[i, :] = line.GetPoint(i)

    t = np.linspace(curv_coor[0], curv_coor[-1], nknots + 2)[1:-1]
    fx = splrep(curv_coor, data[:, 0], k=4, t=t)
    fy = splrep(curv_coor, data[:, 1], k=4, t=t)
    fz = splrep(curv_coor, data[:, 2], k=4, t=t)

    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    data = np.zeros((len(curv_coor), 4))
    data[:, 0] = fx_
    data[:, 1] = fy_
    data[:, 2] = fz_
    data[:, 3] = MISR[:, 0]

    header = ["X", "Y", "Z", radiusArrayName]
    line = data_to_vtkPolyData(data, header)

    # Let vmtk compute curve attributes
    if get_stats:
        line = vmtk_centerline_geometry(line, smooth=False)

    if get_curv:
        # Analytical curvature
        dlsfx = splev(curv_coor, fx, der=1)
        dlsfy = splev(curv_coor, fy, der=1)
        dlsfz = splev(curv_coor, fz, der=1)

        ddlsfx = splev(curv_coor, fx, der=2)
        ddlsfy = splev(curv_coor, fy, der=2)
        ddlsfz = splev(curv_coor, fz, der=2)

        C1xC2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
        C1xC2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
        C1xC2_3 = ddlsfy * dlsfx - ddlsfx * dlsfy

        curvature = np.sqrt(C1xC2_1 ** 2 + C1xC2_2 ** 2 + C1xC2_3 ** 2) / \
                            (dlsfx ** 2 + dlsfy ** 2 + dlsfz ** 2) ** 1.5

        return line, curvature
    else:
        return line


def prepare_surface(model_path, parameters):
    """
    Clean and check connectivity of surface.
    Capps or uncapps surface model at inlet and outlets.

    Args:
        model_path (str): Path to model.
        parameters (dict): Contains surface model information.

    Returns:
        open_surface (vtkPolyData): Open surface model.
    Returns:
        capped_surface (vtkPolyData): Closed surface model.
    """
    # Clean surface
    surface = read_polydata(model_path)
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

    # Check connectivity and only choose the surface with the largest area
    if "check_surface" not in parameters.keys():
        connected_surface = gat_connectivity(surface, mode="Largest")
        if connected_surface.GetNumberOfPoints() != surface.GetNumberOfPoints():
            WritePolyData(surface, model_path.replace(".vtp", "_test.vtp"))

    # Get a capped and uncapped version of the surface
    if is_surface_capped(surface):
        open_surface = uncapp_surface(surface)
        capped_surface = surface
    else:
        open_surface = surface
        capped_surface = capp_surface(surface)

    return open_surface, capped_surface


def prepare_voronoi_diagram(model_path, voronoi_path, voronoi_smoothed_path,
                            smooth, smooth_factor, centerlines):
    """
    Compute and smooth voronoi diagram of surface model.

    Args:
        model_path (str): Path to surface model.
        voronoi_path (str): (Save)path to voronoi diagram.
        voronoi_smoothed_path (str): (Save)path to smoothed voronoi diagram.
        smooth (bool): Voronoi is smoothed if True.
        smooth_factor (float): Smoothing factor for voronoi smoothing.
        centerlines (vtkPolyData): Centerlines throughout geometry.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram of surface.
    """
    # Smooth voronoi diagram
    voronoi = make_voronoi_diagram(model_path, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi_smoothed = smooth_voronoi_diagram(voronoi, centerlines, smooth_factor)
        write_polydata(voronoi_smoothed, voronoi_smoothed_path)
        surface_smoothed = create_new_surface(voronoi_smoothed)
        write_polydata(surface_smoothed, model_smoothed_path)
    else:
        voronoi_smoothed = read_polydata(voronoi_smoothed_path)

    # Use smoothed voronoi or not
    voronoi = voronoi_smoothed if smooth else voronoi

    return voronoi

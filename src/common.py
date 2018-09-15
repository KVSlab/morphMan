import vtk
import numpy as np
import numpy.linalg as la
import sys
import math
import operator

from vmtk import vtkvmtk, vmtkscripts
from vtk.util import numpy_support
from os import path, listdir
from scipy.signal import argrelextrema, gaussian, resample
from scipy.interpolate import splrep, splev
from scipy.ndimage.filters import gaussian_filter as gauss

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
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
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
        input_data (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
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


def get_path_names(input_filepath):
    """Takes the input folderpath as argument, and returns the name of the case name, and
    the path to the parent directory

    Args:
        input_filepath (str): Input filepath

    Returns:
        surface_base_path (str): Name of the parent directory
    """

    surface_name = input_filepath.split(path.sep)[-1].split(".")[0]
    surface_folder_path = path.dirname(input_filepath)
    surface_base_path = path.join(surface_folder_path, surface_name)

    return surface_base_path



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


def smooth_voronoi_diagram(voronoi, centerlines, smoothing_factor,
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
    threshold = get_array(radiusArrayName, centerlines) * (1 - smoothing_factor)

    # Do not smooth inlet and outlets, set threshold to 0
    start = 0
    end = 0
    for i in range(centerlines.GetNumberOfLines()):
        end_ = extract_single_line(centerlines, i).GetNumberOfPoints() - 1
        end += end_
        threshold[start:start+5] = -1
        threshold[end-5:end] = -1
        start += end_ + 1
        end += 1

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

        if distance(point, cl_point) > 2 * threshold[id_] / (1 - smoothing_factor):
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
 #                      centerlines=None):
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
    duplicate points, and/or
    removing unused points
    and/or removing degenerate cells.

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
    """Triangulate surface"""
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
    edges = get_feature_edges(surface)

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
    points = np.zeros((region_array.shape[0], 3))
    for i in range(region_array.shape[0]):
        points[i,] = outputs.GetPoint(i)

    # Get area and center
    area = []
    center = []
    for i in range(int(region_array.max()) + 1):
        # Extract points for this opening
        tmp_points = points[region_array[:, 0] == i, :]

        # Extract the surface edge
        boundary_points = threshold(outputs, "RegionId", lower=i-0.1, upper=i+0.1,
                                    type="between", source=0)

        # Create surface for computing area of opening
        delaunay = vtk.vtkDelaunay2D()
        delaunay.SetInputData(boundary_points)
        delaunay.Update()

        # Add quanteties
        area.append(compute_area(delaunay.GetOutput()))
        center.append(np.mean(tmp_points, axis=0))

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


def remove_distant_points(voronoi, centerline):
    N = voronoi.GetNumberOfPoints()
    newVoronoi = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    radius = np.zeros(N)

    locator = get_locator(centerline)
    get_data = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1
    limit = get_data(0)
    limit = limit * 10

    count = 0
    for i in range(N):
        point = voronoi.GetPoint(i)
        ID = locator.FindClosestPoint(point)
        cl_point = centerline.GetPoint(ID)
        dist = distance(point, cl_point)
        if dist/3 > get_data(i) or get_data(i) > limit:
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


def compute_centerlines(inlet, outlet, filepath, surface, resampling=1,
                        smooth=False, num_iter=100, smooth_factor=0.1,
                        endPoint=1, method="pointlist", recompute=False):
    if path.isfile(str(filepath)) and not recompute:  # Filepath might be None
        return read_polydata(filepath)

    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.SeedSelectorName = method
    centerlines.AppendEndPoints = endPoint
    centerlines.Resampling = 1
    centerlines.ResamplingStepLength = resampling
    centerlines.SourcePoints = inlet
    centerlines.TargetPoints = outlet
    centerlines.Execute()
    centerlines = centerlines.Centerlines

    if smooth:
        centerlineSmoothing = vmtkscripts.vmtkCenterlineSmoothing()
        centerlineSmoothing.SetInputData(centerlines)
        centerlineSmoothing.SetNumberOfSmoothingIterations(num_iter)
        centerlineSmoothing.SetSmoothingFactor(smooth_factor)
        centerlineSmoothing.Update()

        centerlines = centerlineSmoothing.GetOutput()

    # Save the computed centerline.
    if filepath is not None:
        write_polydata(centerlines, filepath)

    return centerlines, centerlines.VoronoiDiagram, centerlines.PoleIds


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
    if not path.isfile(folder + "_info.txt"):
        return {}

    # Get text
    f = open(folder + "_info.txt", "r")
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
    f = open(folder + "info.txt", "w")
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


#def make_centerline(ifile, ofile, length=0.1, it=100, factor=0.1, in_out=None,
#                   smooth=False, resampling=False, newpoints=None,
#                   recompute=False, store_points=False, endpoints=0):
#    """
#    A general centerline command. If a centerline file with the same file
#    name alread exists, then the file is just read. To overwrite this set
#    recompute to True. If recomputed is to True then it uses the exsisting points
#    from the old centerline file and no interaction with the interface is
#    needed. Further on one can choose witch points you want to include by
#    giving in_out. The first element in the list is the source point, and if -1
#    is given, is chooses the old inlet, else it chooses outletX, where X is the
#    outlet ID.
#
#    Args:
#        ifile (vtkPolyData): Input surface model.
#        ofile (str): Name of output file.
#        length (float): Number in (0,1], the resampling step.
#        it(int):  Number of iterations of smoothing
#        factor (float): The smoothening factor
#        smooth (bool): Used to turn on/off smoothening
#        resampling (bool): Used to turn in/off resampling
#        in_out (list): Outlet/inlet points that should be used when recompute
#        newpoints (list): Outlet/inlet points that should be used when recompute
#    """
#
#    # Check if ifile exsists
#    if not path.exists(ifile):
#        print("The input file: %s does not exsists!" % ifile)
#        sys.exit(0)
#
#    # Check if it already exsists or if it is to be recomputed
#    if not path.exists(ofile) or recompute:
#        # If recomputed use the old source and target points
#        parameters = get_parameters(path.dirname(ifile))
#
#        vmtkcenterlines = vmtkscripts.vmtkCenterlines()
#        vmtkcenterlines.Surface = read_polydata(ifile)
#
#        if newpoints is not None:
#            inlet = newpoints[-1]
#            points_ = newpoints[:-1]
#
#        elif "inlet" in parameters:
#            if in_out is None:
#                inlet = parameters["inlet"]
#                out = [k for k in list(parameters.keys()) if "outlet" in k and len(k) < 12]
#                out.sort()
#                points_ = [parameters[p] for p in out]
#            else:
#                inlet = parameters["inlet"] if in_out[0] == -1 else parameters["outlet%s" % in_out[0]]
#                points_ = [parameters["outlet%s" % i] for i in in_out[1:]]
#
#        if newpoints is not None or "inlet" in parameters:
#            # Extract outlet points
#            outlets = []
#            for p in points_:
#                for p_ in p:
#                    outlets.append(p_)
#
#            vmtkcenterlines.SeedSelectorName = "pointlist"
#            vmtkcenterlines.SourcePoints = inlet
#            vmtkcenterlines.TargetPoints = outlets
#        else:
#            vmtkcenterlines.SeedSelectorName = "pickpoint"
#
#        # Add resampling
#        if resampling:
#            vmtkcenterlines.Resampling = 1
#            vmtkcenterlines.ResamplingStepLength = length
#
#        # Exectue command and save centerline
#        vmtkcenterlines.AppendEndPoints = endpoints
#        vmtkcenterlines.Execute()
#
#        centerline = vmtkcenterlines.Centerlines
#        write_polydata(centerline, ofile)
#
#        # Add smoothing
#        if smooth:
#            centerlineSmoothing = vmtkscripts.vmtkCenterlineSmoothing()
#            centerlineSmoothing.Centerlines = centerline
#            centerlineSmoothing.NumberOfSmoothingIterations = it
#            centerlineSmoothing.SmoothingFactor = factor
#            centerlineSmoothing.Execute()
#            centerline = centerlinesSmooth.Centerlines
#
#        # If the points are not already stored, do it now
#        if store_points:
#            for i in range(centerline.GetNumberOfLines()):
#                tmp_line = extract_single_line(centerline, i)
#                tmp_N = tmp_line.GetNumberOfPoints()
#                parameters["outlet%s" % i] = tmp_line.GetPoint(tmp_N - 1)
#            parameters["inlet"] = tmp_line.GetPoint(0)
#            write_parameters(parameters, basedir)
#    else:
#        centerline = read_polydata(ofile)
#
#    return centerline


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
    geometry.FernetTangentArrayName = "FernetTangent"
    geometry.FernetNormalArrayName = "FernetNormal"
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


def spline_centerline(line, get_curv=False, isline=False,
                    nknots=50, get_stats=True, get_misr=True):
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
        get_misr (bool): Determines if MISR values are computed or not.

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
    if get_misr:
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

    if get_misr:
        data = np.zeros((len(curv_coor), 4))
        data[:, 3] = MISR[:, 0]
        header = ["X", "Y", "Z", radiusArrayName]
    else:
        data = np.zeros((len(curv_coor), 3))
        header = ["X", "Y", "Z"]

    data[:, 0] = fx_
    data[:, 1] = fy_
    data[:, 2] = fz_

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


def prepare_surface(surface_path, surface_capped_path, parameters):
    """
    Clean and check connectivity of surface.
    Capps or uncapps surface at inlet and outlets.

    Args:
        surface_path (str): Path to surface.
        parameters (dict): Contains surface information.

    Returns:
        open_surface (vtkPolyData): Open surface.
    Returns:
        capped_surface (vtkPolyData): Closed surface.
    """
    # Check if surface path exists
    if not path.exists(surface_path):
        RuntimeError("Could not find the file: {}".format(surface_path))

    # Clean surface
    surface = read_polydata(surface_path)
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

    # Check connectivity and only choose the surface with the largest area
    if "check_surface" not in parameters.keys():
        connected_surface = get_connectivity(surface, mode="Largest")
        if connected_surface.GetNumberOfPoints() != surface.GetNumberOfPoints():
            WritePolyData(surface, surface_path.replace(".vtp", "_unconnected.vtp"))
            WritePolyData(connected_surface, surface_path)
            surface = connected_surface

        parameters["check_surface"] = True
        write_parameters(parameters, path.dirname(surface_path))

    # Get a capped and uncapped version of the surface
    if is_surface_capped(surface):
        open_surface = uncapp_surface(surface)
        capped_surface = surface
        write_polydata(capped_surface, surface_capped_path)
        write_polydata(open_surface, surface_path)
    else:
        open_surface = surface
        if path.exists(surface_capped_path):
            capped_surface = read_polydata(surface_capped_path)
        else:
            capped_surface = capp_surface(surface)
            write_polydata(capped_surface, surface_capped_path)

    return open_surface, capped_surface


def prepare_voronoi_diagram(surface, capped_surface, centerlines, base_path,
                            smooth, smooth_factor, no_smooth, no_smooth_point):
    """
    Compute and smooth voronoi diagram of surface model.

    Args:
        surface (polydata): Surface model to create a Voronoi diagram of.
        voronoi_path (str): (Save)path to voronoi diagram.
        voronoi_smoothed_path (str): (Save)path to smoothed voronoi diagram.
        smooth (bool): Voronoi is smoothed if True.
        smooth_factor (float): Smoothing factor for voronoi smoothing.
        centerlines (vtkPolyData): Centerlines throughout geometry.
        no_smooth_cl (vtkPolyData): Centerline of region which should not be smoothed.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram of surface.
    """
    # Check if a region should not be smoothed
    if smooth and no_smooth:
        no_smooth_centerline_path = base_path + "_centerlines_no_smooth.vtp"
        # Get inlet and outlets
        inlet = centerlines.GetPoint(0)
        outlets = []
        if no_smooth_point is None:
            seed_selector = vmtkPickPointSeedSelector()
            seed_selector.SetSurface(capped_surface)
            seed_selector.Mode = "no_smooth"
            seed_selector.Execute()
            point_ids = SeedSelector.GetTargetSeedIds()
            outlets = [surface.GetPoint(point_ids.GetId(i)) for i in range(point_ids.GetNumberOfIds())]
        else:
            locator = get_locator(surface)
            for i in range(len(no_smooth_point) // 3):
                outlets.append(surface.GetPoint(locator.FindClosestPoint(no_smooth_point[3*i:3*(i+1)])))

        no_smooth_centerlines = compute_centerlines(inlet, outlets,
                                                    no_smooth_centerline_path,
                                                    capped_surface, resampling=0.1,
                                                    smooth=False)
        no_smooth_segments = []
        for i in range(no_smooth_centerlines.GetNumberOfLines()):
            tmp_line = extract_single_line(no_smooth_centerlines, i)
            div_ids = []
            for j in range(centerlines.GetNumberOfLines()):
                div_ids.append(centerline_div(tmp_line, extract_single_line(centerlines, j)))
            div_id = max(div_ids)
            no_smooth_segments.append(extract_single_line(tmp_line, div_id))

        no_smooth_centerlines = merge_data(no_smooth_segments)

    else:
        no_smooth_cl = None

    # Create Voronoi diagram
    voronoi = make_voronoi_diagram(surface, base_path + "_voronoi.vtp")

    voronoi_smoothed_path = base_path + "_voronoi_smoothed.vtp"
    if not path.exists(voronoi_smoothed_path) and smooth:
        # Smooth voronoi diagram
        voronoi = smooth_voronoi_diagram(voronoi, centerlines, smooth_factor, no_smooth_cl)
        write_polydata(voronoi_smoothed, voronoi_smoothed_path)

        # Create new surface from the smoothed Voronoi
        surface_smoothed = create_new_surface(voronoi_smoothed)
        write_polydata(surface_smoothed, surface_smoothed_path)
    elif smooth:
        voronoi = read_polydata(voronoi_smoothed_path)

    return voronoi


def vtk_box(bounds):
    """Returns a vtk box object based on the bounds

    Args:
        bounds (list): Bounds x1, x2, y1, y2, z1, z2

    Returns:
        box (vtkBox): A vtk box
    """
    box = vtk.vtkBox()
    box.SetBounds(bounds)

    return box


def clip_polydata(surface, cutter):
    """Clip the inpute vtkPolyData object with a cutter function (plane, box, etc)

    Args:
        surface (vtkPolyData): Input vtkPolyData for clipping
        cutter (vtkBox, vtkPlane): Function for cutting the polydata
    """
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(surface)
    clipper.SetClipFunction(cutter)
    clipper.GenerateClipScalarsOn()
    clipper.Update()

    return clipper.GetOutput()


def prepare_surface_output(surface, original_surface, new_centerline, output_filepath,
                           test_merge=False, rotated=False, original_centerline=None,
                           margin=1.1, hight_ratio=0.1):
    # Get planes if outlets of the original surface 
    boundary_edges = get_feature_edges(original_surface)
    boundary_connectivity = get_connectivity(boundary_edges)

    vtk_array = boundary_connectivity.GetPointData().GetArray("RegionId")
    vtk_points = boundary_connectivity.GetPoints().GetData()
    region_id = numpy_support.vtk_to_numpy(vtk_array)
    points = numpy_support.vtk_to_numpy(vtk_points)

    # Get information from the original geometry
    for i in range(region_id.max()):
        # Get relevant points
        tmp_points = points[region_id == i]

        # Get normal
        tmp_normal = np.cross(tmp_points[0] - tmp_points[-1],
                              tmp_points[0] - tmp_points[tmp_points.shape[0] // 2])
        normal = tmp_normal / np.sqrt(np.sum(tmp_normal**2))

        # Get Center
        center = np.mean(tmp_points, axis=0)

        # FIXME: How to deal with the rotated branches when clipping the outlets?
        if rotated:
            pass
            # Translate the center and bounds by comparing the endpoints of the centerlines
            # Rotate the normal and bounds compared to the direction of the vector of the last
            # two points in the centerline. (use center as origo)

        tmp_points = tmp_points - center
        vec = np.eye(3)
        vec[:, 0] = normal
        vec[:, 1] = (tmp_points[0] - center)
        vec[:, 2] = tmp_points[tmp_points.shape[0] // 2] - center
        vec[:, 1] = vec[:, 1] / np.sqrt(np.sum(vec[:, 1]**2))
        vec[:, 2] = vec[:, 2] / np.sqrt(np.sum(vec[:, 2]**2))

        # Make normal z-axis
        R = gram_schmidt(vec)
        R_inv = np.linalg.inv(R)

        tmp_points = np.dot((tmp_points), R)

        # Get bounds
        tmp_points_vtk = vtk.vtkPoints()
        [tmp_points_vtk.InsertNextPoint(tmp_points[k]) for k in range(tmp_points.shape[0])]
        bound = list(tmp_points_vtk.GetBounds())

        # TODO: Test direction

        # Create box with safety margin of 10 % and hight of 0.1
        bound[0] = 0
        bound[1] = abs(bound[2] - bound[3])*hight_ratio
        bound[2] *= margin
        bound[3] *= margin
        bound[4] *= margin
        bound[5] *= margin

        # Rotate 
        rotation_axis = np.array([0, normal[2], -normal[1]]) / np.sqrt(normal[1]**2 + normal[2]**2)
        angle = -np.arccos(normal[0]) * 180 / np.pi

        # Create an implicit function for clipping
        box = vtk_box(bound)

        # Clip the surface
        transform = vtk.vtkTransform()
        transform.Translate(-center[0], -center[1], -center[2])
        translationFilter = vtk.vtkTransformPolyDataFilter()
        translationFilter.SetTransform(transform)
        translationFilter.SetInputData(surface)
        translationFilter.Update()
        surface = translationFilter.GetOutput()

        transform = vtk.vtkTransform()
        transform.RotateWXYZ(-angle, rotation_axis)
        translationFilter = vtk.vtkTransformPolyDataFilter()
        translationFilter.SetTransform(transform)
        translationFilter.SetInputData(surface)
        translationFilter.Update()
        surface = translationFilter.GetOutput()

        surface = clip_polydata(surface, box)

        transform = vtk.vtkTransform()
        transform.RotateWXYZ(angle, rotation_axis)
        translationFilter = vtk.vtkTransformPolyDataFilter()
        translationFilter.SetTransform(transform)
        translationFilter.SetInputData(surface)
        translationFilter.Update()
        surface = translationFilter.GetOutput()

        transform = vtk.vtkTransform()
        transform.Translate(center[0], center[1], center[2])
        translationFilter = vtk.vtkTransformPolyDataFilter()
        translationFilter.SetTransform(transform)
        translationFilter.SetInputData(surface)
        translationFilter.Update()
        surface = translationFilter.GetOutput()

        # Check connectivity and only choose the surface with the largest area
        surface = get_connectivity(surface, mode="Largest")

        write_polydata(surface, "test_rotated_surface.vtp")
        # NOTE: Should test for number of outlets after each clip


    # Perform a 'light' smoothing to obtain a nicer surface
    surface = vmtk_surface_smoother(surface, method="laplace", iterations=100)

    # Clean surface
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

    # Capped surface
    capped_surface = capp_surface(surface)

    if test_merge:
        surface = check_if_surface_is_merged(capped_surface, new_centerline,
                                             output_filepath)

    return surface


#def sort_centerlines(centerlines_complete):
#    """
#    Sort the complete set of centerlines
#    by placing the longest centerline first.
#
#    Args:
#        centerlines_complete (vtkPolyData): Complete set of centerlines.
#
#    Returns:
#        centerlines_in_order (vtkPolyData): Comlete set of sorted centerlines.
#    """
#    lines = []
#    n = centerlines_complete.GetNumberOfCells()
#    for i in range(n):
#        lines.append(extract_single_line(centerlines_complete, i))
#
#    longest = [lines[0]]
#    lenlong = get_curvilinear_coordinate(longest[0])
#    for i in range(1,n):
#        tmplong = get_curvilinear_coordinate(lines[i])
#        if len(tmplong) > len(lenlong):
#            lenlong = tmplong
#            longest.insert(0, lines[i])
#        else:
#            longest.append(lines[i])
#
#    centerlines_in_order = merge_data(longest)
#
#    return centerlines_in_order


def check_if_surface_is_merged(surface, centerlines, output_filepath):
    """
    Clean and check surface for overlapping regions.

    Args:
        surface (vtkPolyData): Surface model.
        centerlines (vtkPolyData): New centerlines.
    """
    # Check if the manipulated centerline and the centerline from the new surface
    # significantly differ, if so it is likely that part of the surface is now merged
    centerlines = vmtk_centerline_resampling(centerlines, length=0.1)
    inlet = centerlines.GetPoint(0)
    outlets = []
    lines_to_compare = []
    for i in range(centerlines.GetNumberOfLines()):
        lines_to_compare.append(extract_single_line(centerlines, i))
        outlets += lines_to_compare[-1].GetPoint(lines_to_compare[-1].GetNumberOfPoints() - 1)

    lines_to_check = compute_centerlines(inlet, outlets, None, surface,
                                               resampling=0.1, recompute=True)
    for i in range(centerlines.GetNumberOfLines()):
        line_to_compare = lines_to_compare[i]
        line_to_check = extract_single_line(lines_to_check, i)

        # Compare distance between points along both centerliens
        N = min([line_to_check.GetNumberOfPoints(), line_to_compare.GetNumberOfPoints()])
        tolerance = get_tolerance(line_to_compare) * 100
        for j in range(N):
            p1 = np.asarray(line_to_check.GetPoint(j))
            p2 = np.asarray(line_to_compare.GetPoint(j))
            dist = distance(p1, p2)
            if dist > tolerance:
                tmp_path = output_filepath.replace(".vtp", "_ERROR_MERGED.vtp")
                write_polydata(surface, tmp_path)
                raise RuntimeError(("\nERROR: Model has most likely overlapping regions. Please check" + \
                            " the surface model {} and provide other parameters for" + \
                            " the manipulation or poly_ball_size.").format(tmp_path))


def get_clipping_points(dirpath, filename):
    """
    Read clipping points from file

    Args:
        dirpath (str): Location of directory.
        filename (str): Name of clipping point file.

    Returns:
        clipping_points (ndarray): Clipping points.
    """
    particles = dirpath + filename
    all_points = np.loadtxt(particles)
    clipping_points = all_points

    return clipping_points


def connect_line(line):
    """
    Create edges between points defining a
    centerline, used to construct a
    discrete line in 3D.

    Args:
        line (vtkPoints): Points defining a centerline.

    Returns:
        line (vtkPolyData): Discrete centerline data.
    """
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
    return line


def move_line_horizontally(patch_cl, ID1, ID2, dx_p1, clip=False, eye=False, side=None):
    """
    Iterate through centerline and move line based on a profile
    for horizontal movement. Includes special treatment of
    opthalmic artery if present.

    Args:
        patch_cl (vtkPolyData): Centerline data.
        ID1 (int): Index of first clipping point.
        ID2 (int): Index of second clipping point.
        dx_p1 (ndarray): Direction to move upstream.
        clip (bool): Determines which part of geometry is being moved, True if siphon.
        eye (bool): Determines presence of opthamlic artery.
        side (str): Determines location relative to the middle of the siphon.

    Returns:
        newline (vtkPolyData): Manipulated centerline.
    """

    centerline_loc = get_locator(patch_cl)
    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if clip:
        if eye:
            l1 = extract_single_line(patch_cl, 0)
            l2 = extract_single_line(patch_cl, 1)
            l3 = extract_single_line(patch_cl, 2)
            test_cl = merge_data([l1, l2, l3])

            ID1 = 0
            ID2 = len(get_curvilinear_coordinate(test_cl))
            idmid = int((ID1 + ID2) / 2.)
        else:
            ID1 = 0
            ID2 = len(get_curvilinear_coordinate(patch_cl))
            idmid = int((ID1 + ID2) / 2.)

        for p in range(patch_cl.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(patch_cl.GetPoint(p))

            if cl_id < idmid:
                dist = dx_p1 * (idmid ** 2 - cl_id ** 2) / (idmid ** 2 - ID1 ** 2)
            else:
                if cl_id <= (ID2 - 1):
                    dist = -dx_p1 * (cl_id - idmid) ** (0.5) / (ID2 - idmid) ** (0.5)
                else:
                    locator = get_locator(test_cl)
                    pp = patch_cl.GetPoint(cl_id)
                    id_main = locator.FindClosestPoint(pp)
                    dist = -dx_p1 * (id_main - idmid) ** (0.5) / (id_main - idmid) ** (0.5)

            patch_point = np.asarray(patch_cl.GetPoint(p))

            if la.norm(patch_point) > 0.1:
                points.InsertNextPoint(patch_point + dist)
                verts.InsertNextCell(1)
                verts.InsertCellPoint(p)

    else:
        for p in range(patch_cl.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(patch_cl.GetPoint(p))

            if side == "right":
                dist = -dx_p1
            elif side == "left":
                dist = dx_p1

            points.InsertNextPoint(np.asarray(patch_cl.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(patch_cl.GetPointData().GetArray(radiusArrayName))

    return newline


def move_points_vertically(line, dx):
    """
    Iterate through centerline points and move line based on a profile
    for vertical movement.

    Args:
        line (vtkPolyData): Centerline data.
        dx (ndarray): Direction to move centerline.

    Returns:
        newline (vtkPolyData): Manipulated centerline.
    """

    centerline_loc = get_locator(line)

    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    ID1 = 0
    ID2 = len(get_curvilinear_coordinate(line))
    for p in range(line.GetNumberOfPoints()):
        cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))

        dist = 4 * dx * (cl_id - ID1) * (ID2 - cl_id) / (ID2 - ID1) ** 2

        points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
        verts.InsertNextCell(1)
        verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(line.GetPointData().GetArray(radiusArrayName))
    return newline


def move_line_vertically(line, dx, ID1_0, clip_ID=None, eye=False):
    """
    Iterate through centerline and move line based on a profile
    for vertical movement. Includes special treatment of
    opthalmic artery if present.

    Args:
        line (vtkPolyData): Centerline data.
        dx (ndarray): Direction to move vertically.
        ID1_0 (int): Index of first clipping point.
        clip_ID (int): Index of opthalmic artery entrance if present.
        eye (bool): Determines presence of opthamlic artery.

    Returns:
        newline (vtkPolyData): Manipulated centerline.
    """
    centerline_loc = get_locator(line)

    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if eye:
        l1 = extract_single_line(line, 0)
        l2 = extract_single_line(line, 1)
        l3 = extract_single_line(line, 2)
        test_cl = merge_data([l1, l2, l3])

        ID1 = I1 = 0
        ID2 = len(get_curvilinear_coordinate(test_cl))
        IDmid = int((ID1 + ID2) / 2.)
        I2 = ID2 - 1

        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))

            if cl_id <= I2:
                dist = 4 * dx * (cl_id - I1) * (I2 - cl_id) / (I2 - I1) ** 2
            else:
                cl_id = clip_ID - ID1_0 + int((ID2 - (clip_ID - ID1_0)) * 0.4)
                dist = 4 * dx * (cl_id - ID1) * (ID2 - cl_id) / (ID2 - ID1) ** 2

            points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)
    else:
        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(line))

        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))
            dist = 4 * dx * (cl_id - ID1) * (ID2 - cl_id) / (ID2 - ID1) ** 2

            points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(line.GetPointData().GetArray(radiusArrayName))

    return newline


def move_perp(n, P, Z, alpha):
    """
    Find directions for manipulation
    in the vertical direction.

    Args:
        n (ndarray): Normal vector to plane through clipping points.
        P (ndarray): Clipping points.
        Z (ndarray): Points along the centerline.
        alpha (float): Extension / Compression factor.

    Returns:
        dZ (ndarray): Directions to points along the  centerline.
        dx (ndarray): Direction to move the centerline.
    """

    p1 = P[0]
    p2 = P[1]

    # Find midpoint and point furthest away
    dist = []
    for z in Z:
        d = la.norm(np.cross((z - p1), (z - p2))) / la.norm(p2 - p1)
        dist.append(d)

    D_id = dist.index(max(dist))
    D = max(dist)
    z_m = Z[D_id]

    # Vector from line to Z_max and projection onto plane
    v = (z_m - p2) - (z_m - p2).dot(p2 - p1) * (p2 - p1) / la.norm(p2 - p1) ** 2
    PV = v - v.dot(n) * n

    # Find distances
    P1 = (z_m - p1).dot(p2 - p1) * (p2 - p1) / la.norm(p2 - p1) ** 2
    P1 = p1 + P1
    V = P1 + v
    PV1 = P1 + PV

    # Move points
    dZ = []
    for i in range(len(Z)):
        dz = np.array(PV) * dist[i] / D * alpha
        dZ.append(Z[i] + dz)

    dx = (PV1 - P1) * alpha

    return dZ, dx


def move_para(n, P, Z, beta):
    """
    Find directions for manipulation
    in the horizontal direction.

    Args:
        n (ndarray): Normal vector to plane through clipping points.
        P (ndarray): Clipping points.
        Z (ndarray): Points along the centerline.
        beta (float): Extension / Compression factor.

    Returns:
        dZ (ndarray): Directions to points along the  centerline.
        zp_min (ndarray): Translation direction in upstream direction.
        zm_min (ndarray): Translation direction in downstream direction.
    """
    p1 = P[0]
    p2 = P[1]

    # Find normal from q
    q = [(p1 + p2) / 2.]
    qp2 = np.array(p2 - q)[0]
    qp1 = np.array(p1 - q)[0]
    s = q[0] - np.cross(qp2, n) * 3

    # Split points based on orientation
    # to q normal
    Z_p = []
    Z_p_dist = []
    Z_m = []
    Z_m_dist = []
    for z in Z:
        d = la.norm(np.cross((z - s), (z - q[0]))) / la.norm(q[0] - s)
        c = np.cross(s - q, z - q)
        if c[0][0] >= 0:
            Z_p.append(z)
            Z_p_dist.append(d)
        else:
            Z_m.append(z)
            Z_m_dist.append(d)

    # Move points
    D = la.norm(p1 - p2) / 2.
    dZ = []
    dZ.append(p1 + qp1 * beta)
    for i in range(len(Z_p)):
        dz = qp1 * Z_p_dist[i] / D * beta
        dZ.append(Z_p[i] + dz)

    for i in range(len(Z_m)):
        dz = qp2 * Z_m_dist[i] / D * beta
        dZ.append(Z_m[i] + dz)
    dZ.append(p2 + qp2 * beta)

    # Check if moved in right direction
    d_0 = la.norm(np.cross(Z_p[0] - q, Z_p[0] - s)) / la.norm(s - q)
    d_1 = la.norm(np.cross(dZ[1] - q, dZ[1] - s)) / la.norm(s - q)
    if d_1 < d_0:
        # Split points based on orientation
        # to q normal
        Z_p = []
        Z_p_dist = []
        Z_m = []
        Z_m_dist = []
        for z in Z:
            d = la.norm(np.cross((z - s), (z - q[0]))) / la.norm(q[0] - s)
            c = -np.cross(s - q, z - q)
            if c[0][0] >= 0:
                Z_p.append(z)
                Z_p_dist.append(d)
            else:
                Z_m.append(z)
                Z_m_dist.append(d)

        # Move points
        D = la.norm(p1 - p2) / 2.
        dZ = []
        dZ.append(p1 + qp1 * beta)
        for i in range(len(Z_p)):
            dz = qp1 * Z_p_dist[i] / D * beta
            dZ.append(Z_p[i] + dz)

        for i in range(len(Z_m)):
            dz = qp2 * Z_m_dist[i] / D * beta
            dZ.append(Z_m[i] + dz)
        dZ.append(p2 + qp2 * beta)

    zpid = Z_p_dist.index(min(Z_p_dist))
    zp_min = Z_p[zpid]
    zmid = Z_m_dist.index(min(Z_m_dist))
    zm_min = Z_m[zmid]

    return dZ, zp_min, zm_min


def best_plane(Z, P):
    """
    Find the least squares plane through
    the points in P and approximating the points
    in Z.

    Args:
        Z (ndarray): Array of points to approximate.
        P (ndarray): Array of points used as constraints.

    Returns:
        n (ndarray): Normal vector to plane.
    """
    # Defined matrices
    b = np.ones(len(Z))
    d = np.ones(len(P))

    # Create complete matrix
    ATA = np.transpose(Z).dot(Z)
    M0 = np.c_[ATA, np.transpose(P)]
    M1 = np.c_[P, np.zeros((len(P), len(P)))]
    M = np.r_[M0, M1]
    Y = np.r_[np.transpose(Z).dot(b), d]

    # Solve system
    x = la.solve(M, Y)
    a = x[0]
    b = x[1]
    c = x[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    # Define plane
    xmin = min(Z, key=operator.itemgetter(1))[0] - 4
    xmax = max(Z, key=operator.itemgetter(1))[0] + 4
    ymin = min(Z, key=operator.itemgetter(1))[1] - 4
    ymax = max(Z, key=operator.itemgetter(1))[1] + 4
    xx, yy = np.meshgrid(np.linspace(xmin, xmax, 15), np.linspace(ymin, ymax, 15))
    zz = (1 - a * xx - b * yy) / float(c)

    return n


def find_closest_point(dx, start, stop, P0, line):
    """
    Find point located closest to a given point P0.
    Searching from start to stop along the centerline.

    Args:
        dx (ndarray): Direction to search for point furthest away.
        start (int): Index to start searching.
        stop (int): Index to stop searching.
        P0 (ndarray): Point to search from.
        line (vtkPolyData): Centerline to search along.

    Returns:
        minP (ndarray): Point located closest to P0.
        minID (int): ID of point located closest to P0.
    """
    a = dx[0]
    b = dx[1]
    c = dx[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    # Define plane
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 100
    xx, yy = np.meshgrid(np.linspace(xmin, xmax, 150), np.linspace(ymin, ymax, 150))
    d = a * P0[0] + b * P0[1] + c * P0[2]
    zz = (d - a * xx - b * yy) / float(c)

    points = []
    for i in range(start, stop):
        p = line.GetPoint(i)
        points.append(np.array(p))

    dist_list = []
    for i, pcl in enumerate(points):
        v = pcl - np.array(P0)
        dist = abs(v.dot(n))
        dist_list.append(dist)

    minID = dist_list.index(min(dist_list)) + start
    minP = points[minID - start]

    return minP, minID


def find_furthest_points(dx, line):
    """
    Find point located furthes away from the line
    spanned of the clipping points p1 and p2.

    Args:
        dx (ndarray): Direction to search for point furthest away.
        line (vtkPolyData): Centerline to search along.

    Returns:
        maxP (ndarray): Point located furthest away.
        maxID (int): ID of point located furthest away.
    """
    P0 = line.GetPoint(0)
    a = dx[0]
    b = dx[1]
    c = dx[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

    # Define plane
    xmin = 0
    xmax = 100
    ymin = 0
    ymax = 100
    xx, yy = np.meshgrid(np.linspace(xmin, xmax, 150), np.linspace(ymin, ymax, 150))
    d = a * P0[0] + b * P0[1] + c * P0[2]
    zz = (d - a * xx - b * yy) / float(c)

    points = []
    for i in range(line.GetNumberOfPoints()):
        p = line.GetPoint(i)
        points.append(np.array(p))

    dist_list = []
    for i, pcl in enumerate(points):
        v = pcl - np.array(P0)
        dist = abs(v.dot(n))
        dist_list.append(dist)

    maxID = dist_list.index(max(dist_list))
    maxP = points[maxID]

    return maxP, maxID


def get_spline_points(line, param, direction, clip_points):
    """
    Pick n uniformly selected points along the
    centerline from point P1 to P2, and move them.

    Args:
        line (vtkPolyData): Longest centerline in geometry.
        param (float): Extension / Compression factor.
        direction (str): Direction to move centerline.
        clip_points (vtkPoints): Clipping points.

    Returns:
        dz (ndarray): Points along the centerline.
        ids (ndarray): IDs of points along centerline.
        dx (ndarray): Direction to move geometry.
    """
    locator = get_locator(line)
    p1 = clip_points.GetPoint(0)
    p2 = clip_points.GetPoint(1)
    ID1 = locator.FindClosestPoint(p1)
    ID2 = locator.FindClosestPoint(p2)
    ID_mid = int((ID1 + ID2) / 2.)
    P = [p1, p2]

    # Select n uniformly spaced points
    n = 10
    points = []
    ids = np.zeros(n)
    dx = 1 / (n + 1.)
    for i in range(1, n + 1):
        ID = int(ID1 + (ID2 - ID1) * i * dx)
        ids[i - 1] = ID
        p = line.GetPoints().GetPoint(ID)
        points.append(np.array([p[0], p[1], p[2]]))

    for i in range(len(P)):
        P[i] = np.array([P[i][0], P[i][1], P[i][2]])

    n = best_plane(points, P)

    if direction == "vertical":
        dz, dx = move_perp(n, P, points, param)
        return dz, ids, dx

    elif direction == "horizont":
        dz, zp, zm = move_para(n, P, points, param)
        return dz, ids


def find_diverging_centerlines(centerlines, end_point):
    """
    Collect centerlines diverging from the longest
    centerline.

    Args:
        centerlines (vtkPolyData): Collection of all centerlines.
        end_point (tuple): End point of relevant centerline.

    Returns:
        div_ids (list): ID where lines diverege.
    Returns:
        div_points (list): Points where lines diverge.
    Returns:
        centerlines (vtkPolyData): Collection of centerlines not diverging.
    Returns:
        div_lines (list): Collection of divering centerlines.
    """
    # Start with longest line
    longest = extract_single_line(centerlines, 0)
    longest = vmtk_centerline_resampling(longest, 0.1)
    longest_locator = get_locator(longest)
    longest_end_id = longest_locator.FindClosestPoint(end_point)

    # Separate lines and divering lines
    lines = [longest]
    div_lines = []
    div_ids = []
    div_points = []
    n = centerlines.GetNumberOfCells() - 1
    tol = 0.40

    # Find diverging lines
    for i in range(n):
        div = False
        line_tmp = extract_single_line(centerlines, i + 1)
        line_tmp = vmtk_centerline_resampling(line_tmp, 0.1)
        stop_id = len(get_curvilinear_coordinate(line_tmp))
        if longest_end_id < stop_id:
            stop_id = longest_end_id
        for j in np.arange(stop_id):
            p_cl = np.asarray(longest.GetPoint(j))
            p_tmp = np.asarray(line_tmp.GetPoint(j))
            dist = distance(p_cl, p_tmp)
            if dist > tol and j < (longest_end_id-20):
                div_lines.append(line_tmp)
                div_ids.append(j)
                div_points.append(p_tmp)
                div = True
                break
        if not div:
            lines.append(line_tmp)

    centerlines = merge_data(lines)

    return div_ids, div_points, centerlines, div_lines


def clip_eyeline(eyeline, clip_start_point, clip_end_ID):
    """
    Clip the opthamlic artery if present.

    Args:
        eyeline (vtkPolyData): Line representing the opthalmic artery centerline.
        clip_star_point (tuple): Point at entrance of opthalmic artery.
        cip_end_ID (int): ID of point at end of opthalmic artery.

    Returns:
        patch_eye (vtkPolyData): Voronoi diagram representing opthalmic artery.
    """

    points = [clip_start_point, eyeline.GetPoint(clip_end_ID)]
    eye_points = vtk.vtkPoints()
    for p in points:
        eye_points.InsertNextPoint(p)

    patch_eye = CreateParentArteryPatches(eyeline, eye_points, siphon=True)

    return patch_eye


def find_ophthalmic_artery(centerlines, clip_pts):
    """
    Method checks if the geometry includes the opthamlic artery.
    Extracts the opthalmic artery if present, and determines its position.

    Args:
        centerlines (vtkPolyData): Complete set of centerlines.
        clip_pts (vtkPoints): Clipping points.

    Returns:
        eye (bool): True if opthamlic artery is present.
        clip_ID (long): ID where opthamlic artery is located along centerline.
        centerlines (vtkPolyData): Complete set of centerlines excluding opthamlic artery.
        eyeline (vtkPolyData): Centerline leading to opthalmic artery.
    """

    # Extract lines:
    lines = []
    n = centerlines.GetNumberOfCells()
    for i in range(n):
        lines.append(extract_single_line(centerlines, i))

    longest = lines[0]
    tol = 0.40

    # Find start and stop IDs along clipped curve
    locator_longest = get_locator(longest)
    p1 = clip_pts[0]
    p2 = clip_pts[1]
    ID1 = locator_longest.FindClosestPoint(p1)
    ID2 = locator_longest.FindClosestPoint(p2)

    eye = False
    index = 1
    for line in lines[1:]:
        locator_checkline = get_locator(line)
        len_check = len(get_curvilinear_coordinate(line))
        if len_check < ID2:
            IDStop = len_check - 1
        else:
            IDStop = ID2

        for i in np.arange(ID1, IDStop):
            p_eye = np.asarray(line.GetPoint(i))
            p_cl = np.asarray(longest.GetPoint(i))
            dist = la.norm(p_eye - p_cl)
            if dist > tol:
                clip_ID = i
                eye = True
                eyeline = lines[index]
                del lines[index]
                centerlines = merge_data(lines)
                break
        if eye:
            break
        index += 1

    if eye:
        return eye, clip_ID, centerlines, eyeline
    else:
        return eye, None, centerlines, None


def get_vtk_clipping_points(line, clipping_points):
    """
    Store clipping points as VTK objects.
    Extract points as tuples and corresponding IDs.

    Args:
        line (vtkPolyData): Line representing longest single centerline.
        clipping_points (ndarray): Array containing two clpping points.

    Returns:
        p1 (tuple): First clipping point.
        p2 (tuple): Second clipping point.
        ID1 (long): ID of first clipping point.
        ID2 (long): ID of second clipping point.
        clip_points (vtkPoints): VTK objects containing the clipping points.
        clipping_points (ndarray): Array containing two clipping points.
    """
    locator = get_locator(line)
    ID1 = locator.FindClosestPoint(clipping_points[0])
    ID2 = locator.FindClosestPoint(clipping_points[1])
    if ID1 > ID2:
        clipping_points = clipping_points[::-1]
        ID1, ID2 = ID2, ID1

    # Set clipping points
    div_points = np.asarray(clipping_points)
    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)
    clip_points = points

    p1 = clip_points.GetPoint(0)
    p2 = clip_points.GetPoint(1)

    return p1, p2, ID1, ID2, clip_points, clipping_points


### The following code is adapted from:
### https://github.com/vmtk/vmtk/tree/master/vmtkApps/CerebralAneurysms/ParentVesselReconstruction
### Written by Marina Piccinelli, and distrubuted within vmtk.

def MaskVoronoiDiagram(voronoi, centerlines):
    numberOfCenterlinesPatches = centerlines.GetNumberOfCells()
    numberOfVoronoiPoints = voronoi.GetNumberOfPoints()

    maskArray = vtk.vtkIntArray()
    maskArray.SetNumberOfComponents(1)
    maskArray.SetNumberOfTuples(numberOfVoronoiPoints)
    maskArray.FillComponent(0, 0)

    for i in range(numberOfCenterlinesPatches):
        tangent, center, centerMISR = compute_patch_end_point_parameters(i, centerlines)
        MaskWithPatch(i, tangent, center, centerMISR, maskArray, centerlines, voronoi)

    return maskArray


def compute_patch_end_point_parameters(id, centerlines):
    point0 = [0.0, 0.0, 0.0]
    point1 = [0.0, 0.0, 0.0]
    tan = [0.0, 0.0, 0.0]
    radius0 = -1

    cell = vtk.vtkGenericCell()
    centerlines.GetCell(id, cell)

    if (id == 0):
        point0 = cell.GetPoints().GetPoint(cell.GetNumberOfPoints() - 1)
        point1 = cell.GetPoints().GetPoint(cell.GetNumberOfPoints() - 2)
        radius0 = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(cell.GetPointId(cell.GetNumberOfPoints() - 1))

    else:
        point0 = cell.GetPoints().GetPoint(0)
        point1 = cell.GetPoints().GetPoint(1)
        radius0 = centerlines.GetPointData().GetArray(radiusArrayName).GetTuple1(cell.GetPointId(0))

    tan[0] = point1[0] - point0[0]
    tan[1] = point1[1] - point0[1]
    tan[2] = point1[2] - point0[2]
    vtk.vtkMath.Normalize(tan)

    return tan, point0, radius0


def MaskWithPatch(id, t, c, r, maskArray, centerlines, voronoi):
    patch = extract_single_line(centerlines, id)

    tubeFunction = vtkvmtk.vtkvmtkPolyBallLine()
    tubeFunction.SetInput(patch)
    tubeFunction.SetPolyBallRadiusArrayName(radiusArrayName)

    lastSphere = vtk.vtkSphere()
    lastSphere.SetRadius(r * 1.5)
    lastSphere.SetCenter(c)

    for i in range(voronoi.GetNumberOfPoints()):
        point = [0.0, 0.0, 0.0]
        voronoiVector = [0.0, 0.0, 0.0]

        voronoi.GetPoint(i, point)
        voronoiVector = [point[j] - c[j] for j in range(3)]
        voronoiVectorDot = vtk.vtkMath.Dot(voronoiVector, t)

        tubevalue = tubeFunction.EvaluateFunction(point)
        spherevalue = lastSphere.EvaluateFunction(point)

        if spherevalue < 0.0 and voronoiVectorDot < 0.0:
            continue
        elif (tubevalue <= 0.0):
            maskArray.SetTuple1(i, 1)


def ComputeNumberOfMaskedPoints(dataArray):
    numberOfPoints = 0
    for i in range(dataArray.GetNumberOfTuples()):
        value = dataArray.GetTuple1(i)
        if value == 1:
            numberOfPoints += 1

    return numberOfPoints


def ExtractMaskedVoronoiPoints(voronoi, maskArray):
    numberOfPoints = ComputeNumberOfMaskedPoints(maskArray)

    maskedVoronoi = vtk.vtkPolyData()
    maskedPoints = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPoints)

    count = 0
    for i in range(voronoi.GetNumberOfPoints()):
        point = [0.0, 0.0, 0.0]
        voronoi.GetPoint(i, point)
        pointRadius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        value = maskArray.GetTuple1(i)
        if (value == 1):
            maskedPoints.InsertNextPoint(point)
            radiusArray.SetTuple1(count, pointRadius)
            cellArray.InsertNextCell(1)
            cellArray.InsertCellPoint(count)
            count += 1

    maskedVoronoi.SetPoints(maskedPoints)
    maskedVoronoi.SetVerts(cellArray)
    maskedVoronoi.GetPointData().AddArray(radiusArray)

    return maskedVoronoi


def CreateParentArteryPatches(parentCenterlines, clipPoints, siphon=False, bif=False):
    numberOfDaughterPatches = parentCenterlines.GetNumberOfCells()
    if siphon:
        clipIds, numberOfPatchedCenterlinesPoints = ExtractPatchesIdsSiphon(parentCenterlines, clipPoints)
    else:
        clipIds, numberOfPatchedCenterlinesPoints = ExtractPatchesIds(parentCenterlines, clipPoints)

    pnt = []

    patchedCenterlines = vtk.vtkPolyData()
    patchedCenterlinesPoints = vtk.vtkPoints()
    patchedCenterlinesCellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPatchedCenterlinesPoints)

    if bif:
        clipIds = sorted(clipIds)

    numberOfCommonPatch = clipIds[0] + 1

    patchedCenterlinesCellArray.InsertNextCell(numberOfCommonPatch)

    count = 0
    line = extract_single_line(parentCenterlines, 0)
    getData = line.GetPointData().GetArray(radiusArrayName).GetTuple1
    for i in range(0, numberOfCommonPatch):
        patchedCenterlinesPoints.InsertNextPoint(line.GetPoint(i))
        patchedCenterlinesCellArray.InsertCellPoint(i)
        radiusArray.SetTuple1(i, getData(i))
        count += 1

    for j in range(numberOfDaughterPatches):
        cell = extract_single_line(parentCenterlines, j)

        getData = cell.GetPointData().GetArray(radiusArrayName).GetTuple1
        numberOfCellPoints = cell.GetNumberOfPoints()
        startId = clipIds[j + 1]

        patchNumberOfPoints = numberOfCellPoints - startId
        patchedCenterlinesCellArray.InsertNextCell(patchNumberOfPoints)

        for i in range(startId, cell.GetNumberOfPoints()):
            point = cell.GetPoint(i)
            patchedCenterlinesPoints.InsertNextPoint(point)
            patchedCenterlinesCellArray.InsertCellPoint(count)
            radiusArray.SetTuple1(count, getData(i))
            count += 1

    patchedCenterlines.SetPoints(patchedCenterlinesPoints)
    patchedCenterlines.SetLines(patchedCenterlinesCellArray)
    patchedCenterlines.GetPointData().AddArray(radiusArray)

    return patchedCenterlines


def ExtractPatchesIdsSiphon(parentCl, clipPts, clipped=False):
    clipIds = []
    numberOfPoints = 0
    N = clipPts.GetNumberOfPoints()

    upstreamPoint = clipPts.GetPoint(0)
    downstreamPoint = clipPts.GetPoint(1)

    for j in range(parentCl.GetNumberOfCells()):
        cellLine = extract_single_line(parentCl, j)
        locator = get_locator(cellLine)

        upId = locator.FindClosestPoint(upstreamPoint)
        downId = locator.FindClosestPoint(downstreamPoint)

        if j == 0:
            if clipped:
                clipIds.append(upId - 1)
                clipIds.append(downId + 1)
            else:
                clipIds.append(upId)
                clipIds.append(downId)
            numberOfPoints += upId + 1
            numberOfPoints += cellLine.GetNumberOfPoints() - downId
        else:
            if clipped:
                clipIds.append(downId + 1)
            else:
                clipIds.append(downId)
            numberOfPoints += cellLine.GetNumberOfPoints() - downId

    return clipIds, numberOfPoints


def ExtractPatchesIds(parentCl, clipPts):
    distance = vtk.vtkMath.Distance2BetweenPoints
    clipIds = []
    numberOfPoints = 0
    N = clipPts.GetNumberOfPoints()

    if N == 3:
        commonPoint = clipPts.GetPoint(0)
        pnt_1 = clipPts.GetPoint(1)
        pnt_2 = clipPts.GetPoint(2)
    else:
        pnt_1 = clipPts.GetPoint(0)
        pnt_2 = clipPts.GetPoint(1)

    for j in range(parentCl.GetNumberOfCells()):
        cellLine = extract_single_line(parentCl, j)
        locator = get_locator(cellLine)

        if j == 0 and N == 3:
            upstreamId = locator.FindClosestPoint(commonPoint)
            clipIds.append(upstreamId)
            numberOfPoints += upstreamId + 1

        ID1 = locator.FindClosestPoint(pnt_1)
        ID2 = locator.FindClosestPoint(pnt_2)

        distance1 = math.sqrt(distance(pnt_1, cellLine.GetPoints().GetPoint(ID1)))
        distance2 = math.sqrt(distance(pnt_2, cellLine.GetPoints().GetPoint(ID2)))

        if distance1 > 1 and distance2 > 1:
            ID = 0
        else:
            ID = ID1 if distance1 < distance2 else ID2

        if N == 2:
            clipIds = [ID1, ID2]
            numberOfPoints = cellLine.GetNumberOfPoints()
        else:
            clipIds.append(ID)
            numberOfPoints += cellLine.GetNumberOfPoints() - ID

    return clipIds, numberOfPoints


def InterpolatePatchCenterlines(patchCenterlines, parentCenterlines,
                                additionalPoint, lower, version):

    if additionalPoint is not None:
        additionalPointIds = []
        for i in range(parentCenterlines.GetNumberOfCells()):
            line = extract_single_line(parentCenterlines, i)
            additionalPointIds.append(line.FindPoint(additionalPoint))
    else:
        additionalPointIds = ["" for i in range(parentCenterlines.GetNumberOfCells())]

    interpolatedLines = vtk.vtkPolyData()
    interpolatedPoints = vtk.vtkPoints()
    interpolatedCellArray = vtk.vtkCellArray()

    pointsInserted = 0
    interpolatedCellArray.Initialize()

    for i in range(parentCenterlines.GetNumberOfCells()):
        startingCell = vtk.vtkGenericCell()
        endingCell = vtk.vtkGenericCell()

        numberOfInterpolationPoints = parentCenterlines.GetCell(i).GetNumberOfPoints()

        patchCenterlines.GetCell(0, startingCell)
        patchCenterlines.GetCell(i + 1, endingCell)

        if version:
            splinePoints = InterpolateSpline(startingCell, endingCell, additionalPoint)
        else:
            splinePoints = InterpolateTwoCells(startingCell, endingCell, \
                                               numberOfInterpolationPoints, \
                                               additionalPointIds[i],
                                               additionalPoint, lower)

        interpolatedCellArray.InsertNextCell(splinePoints.GetNumberOfPoints())
        for j in range(splinePoints.GetNumberOfPoints()):
            interpolatedPoints.InsertNextPoint(splinePoints.GetPoint(j))
            interpolatedCellArray.InsertCellPoint(pointsInserted + j)
        pointsInserted += splinePoints.GetNumberOfPoints()

    interpolatedLines.SetPoints(interpolatedPoints)
    interpolatedLines.SetLines(interpolatedCellArray)

    attributeFilter = vtkvmtk.vtkvmtkCenterlineAttributesFilter()
    attributeFilter.SetInputData(interpolatedLines)
    attributeFilter.SetAbscissasArrayName(abscissasArrayName)
    attributeFilter.SetParallelTransportNormalsArrayName(parallelTransportNormalsArrayName)
    attributeFilter.Update()

    attributeInterpolatedLines = attributeFilter.GetOutput()

    return attributeInterpolatedLines


def InterpolateSpline(startCell, endCell, additionalPoint):
    # If the centerline does not pass the bifurcation, return the centerline
    if startCell.GetPoints().GetPoint(0) == endCell.GetPoints().GetPoint(0):
        return endCell.GetPoints()

    # Get number of cells
    points = []
    num_start = startCell.GetNumberOfPoints()
    num_end = endCell.GetNumberOfPoints()
    get_startCell = startCell.GetPoints()
    get_endCell = endCell.GetPoints()

    points = []
    n = 5
    N = 100
    num_centerline_points = 3

    for i in range(num_centerline_points - 1, -1, -1):
        points.append(get_startCell.GetPoint(num_start - n * i - 1))

    if additionalPoint is not None:
        points.append(additionalPoint)

    for i in range(num_centerline_points):
        points.append(get_endCell.GetPoint(i * n))

    curv_coor = np.zeros(len(points))
    for i in range(len(points) - 1):
        curv_coor[i + 1] = curv_coor[i] + math.sqrt(distance(points[i], points[i + 1]))

    points = np.asarray(points)

    fx = splrep(curv_coor, points[:, 0], k=3)
    fy = splrep(curv_coor, points[:, 1], k=3)
    fz = splrep(curv_coor, points[:, 2], k=3)

    curv_coor = np.linspace(curv_coor[0], curv_coor[-1], N)
    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    tmp = []
    for i in range(num_start - n * num_centerline_points):
        tmp.append(get_startCell.GetPoint(i))

    for j in range(N):
        tmp.append([fx_[j], fy_[j], fz_[j]])

    for k in range(n * num_centerline_points, num_end):
        tmp.append(get_endCell.GetPoint(k))

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(tmp))
    for l in range(len(tmp)):
        points.SetPoint(l, tmp[l])

    return points


def InterpolateTwoCells(startCell, endCell, numberOfSplinePoints, additionalPointId,
                        additionalPoint, type):
    points = vtk.vtkPoints()
    xspline = vtk.vtkCardinalSpline()
    yspline = vtk.vtkCardinalSpline()
    zspline = vtk.vtkCardinalSpline()

    numberOfStartCellPoints = startCell.GetNumberOfPoints()
    numberOfEndCellPoints = endCell.GetNumberOfPoints()
    endCellFirstId = numberOfSplinePoints - numberOfEndCellPoints

    for i in range(numberOfStartCellPoints):
        point = startCell.GetPoints().GetPoint(i)
        xspline.AddPoint(float(i), point[0])
        yspline.AddPoint(float(i), point[1])
        zspline.AddPoint(float(i), point[2])

    if additionalPoint is not None:
        xspline.AddPoint(float(additionalPointId), additionalPoint[0])
        yspline.AddPoint(float(additionalPointId), additionalPoint[1])
        zspline.AddPoint(float(additionalPointId), additionalPoint[2])

    for i in range(numberOfEndCellPoints):
        point = endCell.GetPoints().GetPoint(i)
        index = float(endCellFirstId + i)
        xspline.AddPoint(index, point[0])
        yspline.AddPoint(index, point[1])
        zspline.AddPoint(index, point[2])

    xspline.Compute()
    yspline.Compute()
    zspline.Compute()

    points.SetNumberOfPoints(numberOfSplinePoints)
    for i in range(numberOfSplinePoints):
        points.SetPoint(i, xspline.Evaluate(float(i)), yspline.Evaluate(float(i)), zspline.Evaluate(float(i)))
    return points


def ExtractCylindricInterpolationVoronoiDiagram(cellId, pointId, cylinderRadius, voronoi, centerlines):
    isInside = 0

    if cellId == 0:
        cylinderTop = centerlines.GetPoint(pointId)
        cylinderCenter = centerlines.GetPoint(pointId - interpolationHalfSize)
        cylinderBottom = centerlines.GetPoint(pointId - 2 * interpolationHalfSize)
    else:
        cylinderTop = centerlines.GetPoint(pointId)
        cylinderCenter = centerlines.GetPoint(pointId + interpolationHalfSize)
        cylinderBottom = centerlines.GetPoint(pointId + 2 * interpolationHalfSize)

    interpolationDataset = vtk.vtkPolyData()
    interpolationDatasetPoints = vtk.vtkPoints()
    interpolationDatasetCellArray = vtk.vtkCellArray()

    maskArray = vtk.vtkIntArray()
    maskArray.SetNumberOfComponents(1)
    maskArray.SetNumberOfTuples(voronoi.GetNumberOfPoints())
    maskArray.FillComponent(0, 0)

    for i in range(voronoi.GetNumberOfPoints()):
        point = voronoi.GetPoint(i)
        isInside = IsPointInsideInterpolationCylinder(point, cylinderTop, cylinderCenter, cylinderBottom, cylinderRadius)

        if (isInside == 1):
            maskArray.SetTuple1(i, 1)

    numberOfInterpolationPoints = ComputeNumberOfMaskedPoints(maskArray)

    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfInterpolationPoints)

    count = 0
    for i in range(voronoi.GetNumberOfPoints()):
        value = maskArray.GetTuple1(i)
        if (value == 1):
            interpolationDatasetPoints.InsertNextPoint(voronoi.GetPoint(i))
            interpolationDatasetCellArray.InsertNextCell(1)
            interpolationDatasetCellArray.InsertCellPoint(count)
            radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
            radiusArray.SetTuple1(count, radius)
            count += 1

    interpolationDataset.SetPoints(interpolationDatasetPoints)
    interpolationDataset.SetVerts(interpolationDatasetCellArray)
    interpolationDataset.GetPointData().AddArray(radiusArray)

    return interpolationDataset


def IsPointInsideInterpolationCylinder(x, t, c, b, r):
    halfheigth = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(b, t)) / 2.0

    xc = [x[i] - c[i] for i in range(len(x))]
    tb = [t[i] - b[i] for i in range(len(t))]

    xcnorm = vtk.vtkMath.Norm(xc)
    vtk.vtkMath.Normalize(xc)
    vtk.vtkMath.Normalize(tb)

    alpha = math.acos(vtk.vtkMath.Dot(xc, tb))

    parallelxc = xcnorm * math.cos(alpha)
    perpendicularxc = xcnorm * math.sin(alpha)

    thetamin = math.atan(r / halfheigth)
    thetamax = thetamin + (math.pi - 2 * thetamin)

    inside = 0
    if thetamin <= alpha <= thetamax:
        if abs(perpendicularxc) <= r:
            inside = 1
    else:
        if abs(parallelxc) <= halfheigth:
            inside = 1

    return inside


def ComputeNumberOfMaskedPoints(dataArray):
    numberOfPoints = 0
    for i in range(dataArray.GetNumberOfTuples()):
        value = dataArray.GetTuple1(i)
        if value == 1:
            numberOfPoints += 1

    return numberOfPoints


def VoronoiDiagramInterpolation(interpolationcellid, id0, id1, voronoiDataset0,
                                voronoiDataset1, centerlines, step,
                                clippingPoints):
    cellLine = extract_single_line(centerlines, interpolationcellid)

    startPoint = clippingPoints.GetPoint(id0)
    endPoint = clippingPoints.GetPoint(id1)

    startId = cellLine.FindPoint(startPoint)
    endId = cellLine.FindPoint(endPoint)

    gapStartId = startId + 1 * step
    gapEndId = endId - 1 * step
    arrivalId = gapEndId + 1 * step
    endSavingInterval = gapEndId + 1 * step

    numberOfGapPoints = int(math.fabs(gapEndId - gapStartId)) + 1
    numberOfInterpolationPoints = voronoiDataset0.GetNumberOfPoints()
    numberOfCenterlinesPoints = cellLine.GetNumberOfPoints()
    numberOfAddedPoints = numberOfGapPoints * numberOfInterpolationPoints

    finalNewVoronoiPoints = vtk.vtkPoints()
    cellArray = vtk.vtkCellArray()
    finalRadiusArray = get_vtk_array(radiusArrayName, 1, numberOfAddedPoints)
    count = 0
    for i in range(numberOfInterpolationPoints):
        voronoiPoint = voronoiDataset0.GetPoint(i)
        voronoiPointRadius = voronoiDataset0.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

        centerlinePointLocator = get_locator(cellLine)

        closestPointId = centerlinePointLocator.FindClosestPoint(voronoiPoint)
        closestPoint = cellLine.GetPoint(closestPointId)

        voronoiVector = [0.0, 0.0, 0.0]
        voronoiVector[0] = voronoiPoint[0] - closestPoint[0]
        voronoiVector[1] = voronoiPoint[1] - closestPoint[1]
        voronoiVector[2] = voronoiPoint[2] - closestPoint[2]
        voronoiVectorNorm = vtk.vtkMath.Norm(voronoiVector)
        rotationAngle = ComputeVoronoiVectorToCenterlineAngle(closestPointId, voronoiVector, cellLine)

        PTPoints = vtk.vtkPoints()

        range_step = 1 if closestPointId < arrivalId else -1

        for j in range(closestPointId, arrivalId, range_step):
            localtangent = [0.0, 0.0, 0.0]
            localnormal = [0.0, 0.0, 0.0]
            newVoronoiVector = [0.0, 0.0, 0.0]
            newVoronoiPoint = [0.0, 0.0, 0.0]

            transform = vtk.vtkTransform()
            point0 = [0.0, 0.0, 0.0]
            point0 = cellLine.GetPoint(j)

            if (j < numberOfCenterlinesPoints - 1):
                point1 = [0.0, 0.0, 0.0]
                cellLine.GetPoint(j + 1, point1)
                localtangent[0] += point1[0] - point0[0]
                localtangent[1] += point1[1] - point0[1]
                localtangent[2] += point1[2] - point0[2]
            if (j > 0):
                point2 = [0.0, 0.0, 0.0]
                cellLine.GetPoint(j - 1, point2)
                localtangent[0] += point0[0] - point2[0]
                localtangent[1] += point0[1] - point2[1]
                localtangent[2] += point0[2] - point2[2]

            localnormal = cellLine.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(j)
            localnormaldot = vtk.vtkMath.Dot(localtangent, localnormal)

            localtangent[0] -= localnormaldot * localnormal[0]
            localtangent[1] -= localnormaldot * localnormal[1]
            localtangent[2] -= localnormaldot * localnormal[2]
            vtk.vtkMath.Normalize(localtangent)

            transform.RotateWXYZ(rotationAngle, localtangent)
            transform.TransformNormal(localnormal, newVoronoiVector)
            vtk.vtkMath.Normalize(newVoronoiVector)

            newVoronoiPoint[0] = point0[0] + voronoiVectorNorm * newVoronoiVector[0]
            newVoronoiPoint[1] = point0[1] + voronoiVectorNorm * newVoronoiVector[1]
            newVoronoiPoint[2] = point0[2] + voronoiVectorNorm * newVoronoiVector[2]

            PTPoints.InsertNextPoint(newVoronoiPoint)

        numberOfPTPoints = PTPoints.GetNumberOfPoints()

        lastPTPoint = PTPoints.GetPoint(PTPoints.GetNumberOfPoints() - 1)

        voronoiPointLocator = get_locator(voronoiDataset1)

        arrivalVoronoiPointId = voronoiPointLocator.FindClosestPoint(lastPTPoint)
        arrivalVoronoiPoint = voronoiDataset1.GetPoint(arrivalVoronoiPointId)
        arrivalVoronoiPointRadius = voronoiDataset1.GetPointData().GetArray(radiusArrayName).GetTuple1(arrivalVoronoiPointId)

        arrivalCenterlinePointLocator = get_locator(cellLine)

        arrivalCenterlineClosestPointId = arrivalCenterlinePointLocator.FindClosestPoint(arrivalVoronoiPoint)
        arrivalCenterlineClosestPoint = cellLine.GetPoint(arrivalCenterlineClosestPointId)

        arrivalVoronoiVector = [0.0, 0.0, 0.0]
        arrivalVoronoiVector[0] = arrivalVoronoiPoint[0] - arrivalCenterlineClosestPoint[0]
        arrivalVoronoiVector[1] = arrivalVoronoiPoint[1] - arrivalCenterlineClosestPoint[1]
        arrivalVoronoiVector[2] = arrivalVoronoiPoint[2] - arrivalCenterlineClosestPoint[2]
        arrivalVoronoiVectorNorm = vtk.vtkMath.Norm(arrivalVoronoiVector)
        radiusArray = ComputeSpline(voronoiPointRadius, arrivalVoronoiPointRadius, numberOfPTPoints)
        vectorNormArray = ComputeSpline(voronoiVectorNorm, arrivalVoronoiVectorNorm, numberOfPTPoints)

        pointsToGap = (gapStartId - closestPointId) * step

        pointId = pointsToGap
        for k in range(gapStartId, endSavingInterval, step):
            ptpoint = PTPoints.GetPoint(pointsToGap)
            clpoint = cellLine.GetPoint(k)

            vector = [0.0, 0.0, 0.0]
            vector[0] = ptpoint[0] - clpoint[0]
            vector[1] = ptpoint[1] - clpoint[1]
            vector[2] = ptpoint[2] - clpoint[2]
            vtk.vtkMath.Normalize(vector)

            norm = vectorNormArray.GetTuple1(pointsToGap)

            newvector = [0.0, 0.0, 0.0]
            newvector[0] = norm * vector[0]
            newvector[1] = norm * vector[1]
            newvector[2] = norm * vector[2]

            newpoint = [0.0, 0.0, 0.0]
            newpoint[0] = clpoint[0] + newvector[0]
            newpoint[1] = clpoint[1] + newvector[1]
            newpoint[2] = clpoint[2] + newvector[2]

            finalNewVoronoiPoints.InsertNextPoint(newpoint)
            cellArray.InsertNextCell(1)
            cellArray.InsertCellPoint(count)
            if pointsToGap > 0:
                finalRadiusArray.SetTuple1(count, radiusArray.GetTuple1(pointsToGap))
            pointsToGap += 1
            count += 1

    return finalNewVoronoiPoints, finalRadiusArray


def ComputeVoronoiVectorToCenterlineAngle(pointId, vector, centerline):
    point0 = centerline.GetPoint(pointId)
    point1 = centerline.GetPoint(pointId + 1)
    point2 = centerline.GetPoint(pointId - 1)

    tangent = [0.0, 0.0, 0.0]
    for i in range(3):
        tangent[i] += point1[i] - point0[i]
        tangent[i] += point0[i] - point2[i]

    ptnnormal = centerline.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(pointId)
    alpha = ComputeAngleBetweenVectors(ptnnormal, tangent, vector)

    return alpha


def ComputeAngleBetweenVectors(normal, tangent, vector):
    # Compute the tangent component orthogonal to normal
    otangent = [0.0, 0.0, 0.0]
    normalDot = vtk.vtkMath.Dot(tangent, normal)
    otangent[0] = tangent[0] - normalDot * normal[0]
    otangent[1] = tangent[1] - normalDot * normal[1]
    otangent[2] = tangent[2] - normalDot * normal[2]
    vtk.vtkMath.Normalize(otangent)

    # Compute the vector component orthogonal to otangent, i.e. parallel to normal
    vtk.vtkMath.Normalize(vector)
    ovector = [0.0, 0.0, 0.0]
    vectorDot = vtk.vtkMath.Dot(vector, otangent)
    ovector[0] = vector[0] - vectorDot * otangent[0]
    ovector[1] = vector[1] - vectorDot * otangent[1]
    ovector[2] = vector[2] - vectorDot * otangent[2]
    vtk.vtkMath.Normalize(ovector)

    theta = vtkvmtk.vtkvmtkMath.AngleBetweenNormals(normal, ovector)
    theta = vtk.vtkMath.DegreesFromRadians(theta)

    cross = [0.0, 0.0, 0.0]
    vtk.vtkMath.Cross(ovector, normal, cross)
    tangentDot = vtk.vtkMath.Dot(otangent, cross)

    if (tangentDot < 0.0):
        theta = -1.0 * theta

    angle = -theta

    return angle


def ComputeSpline(startValue, endValue, numberOfPoints):
    averageValue = (startValue + endValue) / 2.0
    averageId = int(numberOfPoints / 2)

    splineArray = vtk.vtkDoubleArray()
    splineArray.SetNumberOfComponents(1)
    splineArray.SetNumberOfTuples(numberOfPoints)
    splineArray.FillComponent(0, 0.0)

    spline = vtk.vtkCardinalSpline()
    spline.AddPoint(float(0), startValue)
    spline.AddPoint(float(averageId), averageValue)
    spline.AddPoint(float(numberOfPoints), endValue)
    spline.Compute()

    for i in range(numberOfPoints):
        splineArray.SetTuple1(i, spline.Evaluate(float(i)))

    return splineArray


def InsertNewVoronoiPoints(oldDataset, newPoints, newArray):
    numberOfDatasetPoints = oldDataset.GetNumberOfPoints()
    numberOfNewPoints = newPoints.GetNumberOfPoints()
    numberOfNewVoronoiPoints = numberOfDatasetPoints + numberOfNewPoints

    newDataset = vtk.vtkPolyData()
    cellArray = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfNewVoronoiPoints)

    for i in range(numberOfDatasetPoints):
        point = [0.0, 0.0, 0.0]
        oldDataset.GetPoint(i, point)
        points.InsertNextPoint(point)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(i)

        value = oldDataset.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        radiusArray.SetTuple1(i, value)

    for i in range(numberOfNewPoints):
        point = [0.0, 0.0, 0.0]
        newPoints.GetPoint(i, point)
        points.InsertNextPoint(point)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(numberOfDatasetPoints + i)

        value = newArray.GetTuple1(i)
        radiusArray.SetTuple1(numberOfDatasetPoints + i, value)

    newDataset.SetPoints(points)
    newDataset.SetVerts(cellArray)
    newDataset.GetPointData().AddArray(radiusArray)

    return newDataset


def interpolate_voronoi_diagram(interpolatedCenterlines, patchCenterlines,
                                clippedVoronoi, clippingPoints, bif, lower,
                                cylinder_factor):
    # Extract clipping points
    clippingPointsArray = clippingPoints[1]
    clippingPoints = clippingPoints[0]

    # Copy the voronoi diagram
    completeVoronoiDiagram = vtk.vtkPolyData()
    completeVoronoiDiagram.DeepCopy(clippedVoronoi)

    numberOfInterpolatedCenterlinesCells = interpolatedCenterlines.GetNumberOfCells()
    for j in range(1, 3):
        interpolationCellId = j - 1
        startId = 0
        endId = j

        startCell = vtk.vtkGenericCell()
        patchCenterlines.GetCell(startId, startCell)

        startCellPointId = startCell.GetPointId(startCell.GetNumberOfPoints() - 1)
        startCellPoint = patchCenterlines.GetPoint(startCellPointId)
        startCellPointRadius = patchCenterlines.GetPointData().GetArray(radiusArrayName)\
                                                        .GetTuple1(startCellPointId)
        startCellPointHalfRadius = startCellPointRadius / cylinder_factor

        startInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId,
                                                    startCellPointId, startCellPointRadius,
                                                    clippedVoronoi, patchCenterlines)
        startHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId,
                                                    startCellPointId, startCellPointHalfRadius,
                                                    clippedVoronoi, patchCenterlines)
        endCell = vtk.vtkGenericCell()
        patchCenterlines.GetCell(endId, endCell)

        endCellPointId = endCell.GetPointId(0)
        endCellPoint = patchCenterlines.GetPoint(endCellPointId)
        endCellPointRadius = patchCenterlines.GetPointData().GetArray(radiusArrayName)\
                                                        .GetTuple1(endCellPointId)
        endCellPointHalfRadius = endCellPointRadius / cylinder_factor
        endInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId, endCellPointId,
                                                    endCellPointRadius, clippedVoronoi,
                                                    patchCenterlines)
        endHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId,
                                                    endCellPointId, endCellPointHalfRadius,
                                                    clippedVoronoi, patchCenterlines)

        # Find and insert new points
        newVoronoiPoints, newVoronoiPointsMISR = VoronoiDiagramInterpolation(interpolationCellId,
                                                    startId, endId, startInterpolationDataset,
                                                    endHalfInterpolationDataset,
                                                    interpolatedCenterlines, 1,
                                                    clippingPoints)
        completeVoronoiDiagram = InsertNewVoronoiPoints(completeVoronoiDiagram, newVoronoiPoints,
                                    newVoronoiPointsMISR)

        newVoronoiPoints, newVoronoiPointsMISR = VoronoiDiagramInterpolation(interpolationCellId,
                                                    endId, startId, endInterpolationDataset,
                                                    startHalfInterpolationDataset,
                                                    interpolatedCenterlines, -1,
                                                    clippingPoints)
        completeVoronoiDiagram = InsertNewVoronoiPoints(completeVoronoiDiagram, newVoronoiPoints,
                                                    newVoronoiPointsMISR)

    if bif is not []:
        for i in range(len(bif) - 1):
            interpolationCellId = 0
            startId_ = 0
            endId_ = 1
            bif_clipped = bif[-1]
            bif_ = bif[i]

            startCell = extract_single_line(bif_clipped, startId_)
            locator = get_locator(startCell)
            startPoint = clippingPoints.GetPoint(1)
            startId = locator.FindClosestPoint(startPoint)
            startPoint = startCell.GetPoint(startId)
            startR = startCell.GetPointData().GetArray(radiusArrayName).GetTuple1(startId)
            startRHalf = startR / cylinder_factor

            endCell = extract_single_line(bif_clipped, endId_)
            locator = get_locator(endCell)
            endPoint = endCell.GetPoint(0)
            endId = locator.FindClosestPoint(endPoint)
            endPoint = endCell.GetPoint(endId)
            endR = endCell.GetPointData().GetArray(radiusArrayName).GetTuple1(endId)
            endRHalf = endR / cylinder_factor

            id1, id2 = get_start_ids(clippingPointsArray, bif_clipped)

            startInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId_,
                                                        startId, startR,
                                                        clippedVoronoi, startCell)
            startHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId_,
                                                        startId, startRHalf,
                                                        clippedVoronoi, startCell)

            endInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId_, endId,
                                                        endR, clippedVoronoi,
                                                        endCell)
            endHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId_, endId,
                                                        endRHalf, clippedVoronoi,
                                                        endCell)
            newVoronoiPoints, newVoronoiPointsMISR = VoronoiDiagramInterpolation(interpolationCellId,
                                                        id1, id2, startInterpolationDataset,
                                                        endHalfInterpolationDataset,
                                                        bif_, 1, clippingPoints)

            completeVoronoiDiagram = InsertNewVoronoiPoints(completeVoronoiDiagram, newVoronoiPoints,
                                                        newVoronoiPointsMISR)
            newVoronoiPoints, newVoronoiPointsMISR = VoronoiDiagramInterpolation(interpolationCellId,
                                                        id2, id1, endInterpolationDataset,
                                                        startHalfInterpolationDataset,
                                                        bif_, -1, clippingPoints)

            completeVoronoiDiagram = InsertNewVoronoiPoints(completeVoronoiDiagram, newVoronoiPoints,
                                                            newVoronoiPointsMISR)

            print("Number of points in Voronoi diagram: %i" % completeVoronoiDiagram.GetNumberOfPoints())

    return completeVoronoiDiagram


def get_start_ids(points, line):
    p1 = points[1]
    p2 = points[2]
    locator = get_locator(line)
    id1 = locator.FindClosestPoint(p1)
    id2 = locator.FindClosestPoint(p2)
    if id1 < id2:
        return 1, 2
    else:
        return 2, 1

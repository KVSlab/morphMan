import vtk
import numpy as np
import numpy.linalg as la
import math
from os import path, makedirs
import sys

from vmtk import vtkvmtk, vmtkscripts
from vtk.util import numpy_support
from scipy.signal import resample
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


def read_polydata(filename, type=None):
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
        reader = vtk.vtkXMinkorporereLStructuredGridReader()
    elif fileType == 'vtr':
        reader = vtk.vtkXMLRectilinearGridReader()
    elif fileType == 'vtu':
        reader = vtk.vtkXMLUnstructuredGridReader()
    elif fileType == "vti":
        reader = vtk.vtkXMLImageDataReader()
    elif fileType == "np" and type == "vtkIdList":
        result = np.load(filename).astype(np.int)
        id_list = vtk.vtkIdList()
        id_list.SetNumberOfIds(result.shape[0])
        for i in range(result.shape[0]):
            id_list.SetId(i, result[i])
        return id_list
    else:
        raise RuntimeError('Unknown file type %s' % fileType)

    # Read
    reader.SetFileName(filename)
    reader.Update()
    polyData = reader.GetOutput()

    return polyData


def write_polydata(input_data, filename, type=None):
    """
    Write the given input data based on the file name extension.

    Args:
        input_data (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Input data.
        filename (str): Save path location.
        type (str): Additional parameter for vtkIdList objects.
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
    elif fileType == "np" and type == "vtkIdList":
        output_data = np.zeros(input_data.GetNumberOfIds())
        for i in range(input_data.GetNumberOfIds()):
            output_data[i] = input_data.GetId(i)
        output_data.dump(filename)
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
        base_path (str): Path to the surface, but without the extension
    """

    surface_name = input_filepath.split(path.sep)[-1].split(".")[0]
    folder = path.dirname(input_filepath)
    base_path = path.join(folder, surface_name)

    return base_path


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
    voronoi.RemoveSubresolutionTetrahedra = 1
    voronoi.Execute()

    write_polydata(voronoi.VoronoiDiagram, filename)

    newVoronoi = voronoi.VoronoiDiagram
    return newVoronoi


def get_tolerance(centerline, n=50):
    """
    Finds tolerance based on
    average length between first N points
    along the input centerline.

    Args:
        centerline (vtkPolyData): Centerline data.
        n (int): Number of points

    Returns:
        tolerance (float): Tolerance value.
    """

    line = extract_single_line(centerline, 0)
    length = get_curvilinear_coordinate(line)
    tolerance = np.mean(length[1:n] - length[:n - 1]) / divergingRatioToSpacingTolerance

    return tolerance


def get_relevant_outlets(surface, base_path):
    """
    Extract relevant outlets of the
    input surface model.

    Args:
        surface(vtkPolyData): Surface model.
        base_path (str): Location of info-file.

    Returns:
        relevant_outlets (list): List of relevant outlet IDs.
    """
    # Check if info exists
    if not path.isfile(base_path + "_info.txt"):
        provide_relevant_outlets(surface, base_path)

    # Open info
    parameters = get_parameters(base_path)
    relevant_outlets = []
    for key, value in list(parameters.items()):
        if key.startswith("relevant_outlet_"):
            relevant_outlets.append(value)

    if relevant_outlets == []:
        relevant_outlets = provide_relevant_outlets(surface, base_path)

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
        smoothing_factor (float): Smoothing factor.
        no_smooth_cl (vktPolyData): Unsmoothed centerline.

    Returns: smoothedDiagram (vtkPolyData): Smoothed voronoi diagram.
    """
    numberOfPoints = voronoi.GetNumberOfPoints()
    threshold = get_array(radiusArrayName, centerlines) * (1 - smoothing_factor)

    # Do not smooth inlet and outlets, set threshold to -1
    start = 0
    end = 0
    for i in range(centerlines.GetNumberOfLines()):
        line = extract_single_line(centerlines, i)
        length = get_curvilinear_coordinate(line)
        end_ = line.GetNumberOfPoints() - 1
        end += end_

        # Point buffer start
        end_id = end_ - np.argmin(np.abs(-(length - length.max()) - threshold[end]))
        start_id = np.argmin(np.abs(length - threshold[start]))

        threshold[start:start + start_id] = -1
        threshold[end - end_id:end] = -1
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
    # w numpy support, n=63000 32.7 ms, wo numpy support 33.4 ms
    # vtk_array = line.GetPointData().GetArray(arrayName)
    # array = numpy_support.vtk_to_numpy(vtk_array)
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
    # vtk_array = line.GetCellData().GetArray(arrayName)
    # array = numpy_support.vtk_to_numpy(vtk_array)
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


def create_new_surface(completeVoronoiDiagram, poly_ball_size=[120, 120, 120]):
    """
    Envelops an input voronoi diagram
    into a new surface model at a
    given resolution determined by
    the poly_ball_size.

    Args:
        completeVoronoiDiagram (vtkPolyData): Voronoi diagram
        poly_ball_size (list): List of dimensional resolution of output model

    Returns:
        envelope (vtkPolyData): Enveloped surface model.
    """
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    modeller.SetInputData(completeVoronoiDiagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(poly_ball_size)
    modeller.Update()

    # Write the new surface
    marchingCube = vtk.vtkMarchingCubes()
    marchingCube.SetInputData(modeller.GetOutput())
    marchingCube.SetValue(0, 0.0)
    marchingCube.Update()
    envelope = marchingCube.GetOutput()

    return envelope


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
    cleaned_surface = clean_surface(surface)
    triangulated_surface = triangulate_surface(cleaned_surface)

    # Select seeds
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(triangulated_surface)
    seed_selector.text = "Please select the two relevant outlets, \'u\' to undo\n"
    seed_selector.Execute()

    aneurysmSeedIds = seed_selector.GetTargetSeedIds()
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
    cleaned_surface = clean_surface(surface)
    triangulated_surface = triangulate_surface(cleaned_surface)

    # Select seeds
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(triangulated_surface)
    seed_selector.text = "Please position the mouse and press space to add the top of the" + \
                         " aneurysm, \'u\' to undo\n"
    seed_selector.Execute()

    aneurysmSeedIds = seed_selector.GetTargetSeedIds()
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

    # Declear dictionary to hold results
    data = {"bif": {}, 0: {}, 1: {}}

    # Find lower clipping point
    N = min(cl1.GetNumberOfPoints(), cl2.GetNumberOfPoints())
    for i in range(0, N):
        point_0 = cl1.GetPoint(i)
        point_1 = cl2.GetPoint(i)
        distance_between_points = distance(point_0, point_1)
        if distance_between_points > tol:
            center = cl1.GetPoint(i)
            r = cl1.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
            break

    end, r_end, id_end = move_past_sphere(cl1, center, r, i, step=-1)
    data["bif"]["end_point"] = end
    data["bif"]["div_point"] = center

    # Find the diverging points for the bifurcation
    # continue further downstream in each direction and stop when
    # a point is closer than tol, then move point MISR * X
    locator = get_locator(centerline_bif)

    for counter, cl in enumerate([cl1, cl2]):
        for i in range(i, cl.GetNumberOfPoints(), 1):
            tmp_point = cl.GetPoint(i)
            closest_point_ID = locator.FindClosestPoint(tmp_point)
            closest_point = centerline_bif.GetPoint(closest_point_ID)
            distance_between_points = distance(tmp_point, closest_point)
            if distance_between_points < tol * 4:
                center = cl.GetPoint(i)
                r = cl1.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
                break

        end, r_end, id_end = move_past_sphere(cl, center, r, i, step=1,
                                              stop=i * 100, X=1)
        data[counter]["end_point"] = end
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


def clean_surface(surface):
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


def get_centers(surface, base_path, flowext=False):
    # Check if info exists
    if flowext or not path.isfile(base_path + "_info.txt"):
        compute_centers(surface, base_path)

    # Open info
    parameters = get_parameters(base_path)
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
        inlet, outlets = compute_centers(surface, base_path)

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
    """
    Compute area of polydata
    """
    mass = vtk.vtkMassProperties()
    mass.SetInputData(surface)

    return mass.GetSurfaceArea()


def clipp_capped_surface(surface, centerlines, clipspheres=0):
    extractor = vmtkscripts.vmtkEndpointExtractor()
    extractor.Centerlines = centerlines
    extractor.RadiusArrayName = radiusArrayName
    extractor.GroupIdsArrayName = groupIDsArrayName
    extractor.BlankingArrayName = branchClippingArrayName
    extractor.NumberOfEndPointSpheres = clipspheres
    extractor.Execute()
    clipped_centerlines = extractor.Centerlines

    clipper = vmtkscripts.vmtkBranchClipper()
    clipper.Surface = surface
    clipper.Centerlines = clipped_centerlines
    clipper.RadiusArrayName = radiusArrayName
    clipper.GroupIdsArrayName = groupIDsArrayName
    clipper.BlankingArrayName = branchClippingArrayName
    clipper.Execute()
    surface = clipper.Surface

    connector = vmtkscripts.vmtkSurfaceConnectivity()
    connector.Surface = surface
    connector.CleanOutput = 1
    connector.Execute()
    surface = connector.Surface

    return surface


def uncapp_surface(surface, gradients_limit=0.15, area_limit=0.3, circleness_limit=3):
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

    # Mark all cells with a gradient magnitude less then gradient_limit
    end_capp_array = gradients_magnitude < gradients_limit
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

    # Only keep outlets with circleness < circleness_limit and area > area_limit
    circleness_IDs = np.where(np.array(circleness) < circleness_limit)
    region_IDs = np.where(np.array(area) > area_limit)
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


def check_if_surface_is_capped(surface):
    # Get boundary cells
    cells = get_feature_edges(surface)
    if cells.GetNumberOfCells() == 0:
        return True, 0
    else:
        outlets = get_connectivity(cells, mode="All")
        number = get_array("RegionId", outlets).max()
        return number == 0, int(number)


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
    if point_radius[argsort[1]] / point_radius[argsort[0]] > 5:
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
            return True, 0
        else:
            return False, cells.GetNumberOfCells()

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
        boundary_points = threshold(outputs, "RegionId", lower=i - 0.1, upper=i + 0.1,
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
        if dist / 3 > get_data(i) or get_data(i) > limit:
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


def compute_centerlines(inlet, outlet, filepath, surface, resampling=1.0, smooth=False,
                        num_iter=100, smooth_factor=0.1, endPoint=1, method="pointlist",
                        recompute=False, voronoi=None, pole_ids=None, base_path=None):
    if path.isfile(str(filepath)) and not recompute:  # Filepath might be None
        if base_path is not None and path.isfile(base_path + "_voronoi.vtp"):
            voronoi = read_polydata(base_path + "_voronoi.vtp")
            pole_ids = read_polydata(base_path + "_pole_ids.np", type="vtkIdList")
        else:
            voronoi = None
            pole_ids = None

        return read_polydata(filepath), voronoi, pole_ids

    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.SeedSelectorName = method
    centerlines.AppendEndPoints = endPoint
    centerlines.Resampling = 1
    centerlines.ResamplingStepLength = resampling
    centerlines.SourcePoints = inlet
    centerlines.TargetPoints = outlet
    if voronoi is not None and pole_ids is not None:
        centerlines.VoronoiDiagram = voronoi
        centerlines.PoleIds = pole_ids
    centerlines.Execute()
    centerlines_output = centerlines.Centerlines

    if smooth:
        centerlineSmoothing = vmtkscripts.vmtkCenterlineSmoothing()
        centerlineSmoothing.SetInputData(centerlines_output)
        centerlineSmoothing.SetNumberOfSmoothingIterations(num_iter)
        centerlineSmoothing.SetSmoothingFactor(smooth_factor)
        centerlineSmoothing.Update()

        centerlines_output = centerlineSmoothing.GetOutput()

    # Save the computed centerline.
    if filepath is not None:
        write_polydata(centerlines_output, filepath)

    voronoi = centerlines.VoronoiDiagram
    pole_ids = centerlines.PoleIds
    if base_path is not None:
        write_polydata(voronoi, base_path + "_voronoi.vtp")
        write_polydata(pole_ids, base_path + "_pole_ids.np", type="vtkIdList")

    return centerlines_output, voronoi, pole_ids


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
    f = open(folder + "_info.txt", "w")
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

    return tempPoint, r, i


def vmtk_surface_smoother(surface, method, iterations=800, passband=1.0, relaxation=0.01):
    smoother = vmtkscripts.vmtkSurfaceSmoothing()
    smoother.Surface = surface
    smoother.NumberOfIterations = iterations
    if method == "laplace":
        smoother.RelaxationFactor = relaxation
    elif method == "taubin":
        smoother.PassBand = passband
    smoother.Method = method
    smoother.Execute()
    surface = smoother.Surface

    return surface


def extract_ica_centerline(base_path, resampling_step):
    input_filepath = base_path + ".vtp"
    ica_centerline_path = base_path + "_ica.vtp"
    centerline_relevant_outlets_path = base_path + "_centerline_relevant_outlets.vtp"
    if path.exists(ica_centerline_path):
        return read_polydata(ica_centerline_path)

    # Prepare surface and identify in/outlets
    surface, capped_surface = prepare_surface(base_path, input_filepath)
    inlet, outlets = get_centers(surface, base_path)
    outlet1, outlet2 = get_relevant_outlets(capped_surface, base_path)
    outlets, outlet1, outlet2 = sort_outlets(outlets, outlet1, outlet2, base_path)

    # Get relevant centerlines
    centerline_relevant_outlets = compute_centerlines(inlet, outlet1 + outlet2,
                                                      centerline_relevant_outlets_path,
                                                      capped_surface,
                                                      resampling=resampling_step)

    # Extract ICA centerline
    tmp_line_1 = extract_single_line(centerline_relevant_outlets, 0)
    tmp_line_2 = extract_single_line(centerline_relevant_outlets, 1)
    line = extract_single_line(centerlines, 0, startID=0, endID=centerline_div(tmp_line_1,
                                                                               tmp_line_2))
    write_polydata(line, ica_centerline_path)
    return line


def sort_outlets(outlets, outlet1, outlet2, dirpath):
    """
    Sort all outlets of the geometry given the two relevant outlets

    Args:
        outlets (list): List of outlet center points.
        outlet1 (list): Point representing first relevant oultet.
        outlet2 (list): Point representing second relevant oultet.
        dirpath (str): Location of info file.

    Returns:
        outlets (list): List of sorted outlet center points.
    Returns:
        outlet1 (list): Point representing first relevant oultet.
    Returns:
        outlet2 (list): Point representing second relevant oultet.
    """
    tmp_outlets = np.array(outlets).reshape(len(outlets) // 3, 3)
    outlet1_index = np.argsort(np.sum((tmp_outlets - outlet1) ** 2, axis=1))[0]
    outlet2_index = np.argsort(np.sum((tmp_outlets - outlet2) ** 2, axis=1))[0]
    tmp_outlets = tmp_outlets.tolist()
    if max(outlet1_index, outlet2_index) == outlet1_index:
        outlet1 = tmp_outlets.pop(outlet1_index)
        outlet2 = tmp_outlets.pop(outlet2_index)
    else:
        outlet2 = tmp_outlets.pop(outlet2_index)
        outlet1 = tmp_outlets.pop(outlet1_index)
    outlet_rest = (np.array(tmp_outlets).flatten()).tolist()
    outlets = outlet1 + outlet2 + outlet_rest
    data = {}
    for i in range(len(outlets) // 3):
        data["outlet" + str(i)] = outlets[3 * i:3 * (i + 1)]
    write_parameters(data, dirpath)

    return outlets, outlet1, outlet2


def vmtk_centerline_resampling(line, length):
    resampler = vmtkscripts.vmtkCenterlineResampling()
    resampler.Centerlines = line
    resampler.Length = length
    resampler.Execute()

    line = resampler.Centerlines

    return line


def str2bool(boolean):
    if boolean.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif boolean.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ValueError('Boolean value expected.')


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

    k2 = (curvature.T * (E1[:, 1] * N[:, 0] - N[:, 1] * E1[:, 0]) /
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


def prepare_surface(base_path, surface_path):
    """
    Clean and check connectivity of surface.
    Capps or uncapps surface at inlet and outlets.

    Args:
        base_path (str): Absolute path to base folder.
        surface_path (str): Path to surface.

    Returns:
        open_surface (vtkPolyData): Open surface.
    Returns:
        capped_surface (vtkPolyData): Closed surface.
    """
    # Check if surface path exists
    surface_capped_path = base_path + "_capped.vtp"
    if not path.exists(surface_path):
        RuntimeError("Could not find the file: {}".format(surface_path))

    # Clean surface
    surface = read_polydata(surface_path)
    surface = clean_surface(surface)
    surface = triangulate_surface(surface)

    # Check connectivity and only choose the surface with the largest area
    parameters = get_parameters(base_path)
    if "check_surface" not in parameters.keys():
        connected_surface = get_connectivity(surface, mode="Largest")
        if connected_surface.GetNumberOfPoints() != surface.GetNumberOfPoints():
            write_polydata(surface, surface_path.replace(".vtp", "_unconnected.vtp"))
            write_polydata(connected_surface, surface_path)
            surface = connected_surface

        parameters["check_surface"] = True
        write_parameters(parameters, base_path)

    # Get a capped and uncapped version of the surface
    cap_bool, num_out = check_if_surface_is_capped(surface)
    if cap_bool:
        open_surface = uncapp_surface(surface)
        cap_bool, num_out = check_if_surface_is_capped(open_surface)
        print(("WARNING: Tried to automagically uncapp the input surface. Uncapped {}" +
               " inlet/outlets in total. If this number if incorrect please provide an" +
               " uncapped surface as input, use the clipp_capped_surface" +
               " method, or vmtksurfaceendclipper.").format(num_out))
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


def prepare_voronoi_diagram(capped_surface, centerlines, base_path,
                            smooth, smooth_factor, no_smooth, no_smooth_point,
                            voronoi, pole_ids):
    """
    Compute and smooth voronoi diagram of surface model.

    Args:
        capped_surface (polydata): Cappedsurface model to create a Voronoi diagram of.
        base_path (str): Absolute path to surface model path.
        voronoi (vtkPolyData): Voronoi diagram.
        pole_ids (vtkIDList): Pole ids of Voronoi diagram.
        smooth (bool): Voronoi is smoothed if True.
        smooth_factor (float): Smoothing factor for voronoi smoothing.
        centerlines (vtkPolyData): Centerlines throughout geometry.
        no_smooth (bool): Part of Voronoi is not smoothed.
        no_smooth_point (vtkPolyData): Point which defines unsmoothed area.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram of surface.
    """
    # Check if a region should not be smoothed
    if smooth and no_smooth:
        no_smooth_path = base_path + "_centerline_no_smooth.vtp"
        if not path.exists(no_smooth_path):
            # Get inlet and outlets
            tol = get_tolerance(centerlines)
            inlet = extract_single_line(centerlines, 0)
            inlet = inlet.GetPoint(0)  # inlet.GetNumberOfPoints() - 1)
            outlets = []
            if no_smooth_point is None:
                seed_selector = vmtkPickPointSeedSelector()
                seed_selector.SetSurface(capped_surface)
                seed_selector.text = "Please place a point on the segments you do not want" + \
                                     " to smooth, e.g. an aneurysm, \'u\' to undo\n"
                seed_selector.Execute()
                point_ids = seed_selector.GetTargetSeedIds()
                for i in range(point_ids.GetNumberOfIds()):
                    outlets += capped_surface.GetPoint(point_ids.GetId(i))
            else:
                locator = get_locator(capped_surface)
                for i in range(len(no_smooth_point) // 3):
                    tmp_id = locator.FindClosestPoint(no_smooth_point[3 * i:3 * (i + 1)])
                    outlets.append(capped_surface.GetPoint(tmp_id))

            # Create the centerline
            no_smooth_centerlines, _, _ = compute_centerlines(inlet, outlets,
                                                              None, capped_surface,
                                                              resampling=0.1,
                                                              smooth=False, voronoi=voronoi,
                                                              pole_ids=pole_ids)

            # Extract the centerline region which diverges from the existing centerlines
            no_smooth_segments = []
            for i in range(no_smooth_centerlines.GetNumberOfLines()):
                tmp_line = extract_single_line(no_smooth_centerlines, i)
                div_ids = []
                for j in range(centerlines.GetNumberOfLines()):
                    div_ids.append(centerline_div(tmp_line, extract_single_line(centerlines, j), tol))
                div_id = max(div_ids)
                no_smooth_segments.append(extract_single_line(tmp_line, 0, startID=div_id))

            no_smooth_cl = merge_data(no_smooth_segments)
            write_polydata(no_smooth_cl, no_smooth_path)
        else:
            no_smooth_cl = read_polydata(no_smooth_path)

    else:
        no_smooth_cl = None

    if voronoi is None:
        voronoi = make_voronoi_diagram(capped_surface, base_path + "_voronoi.vtp")

    # Smooth voronoi
    voronoi_smoothed_path = base_path + "_voronoi_smoothed.vtp"
    surface_smoothed_path = base_path + "_smoothed.vtp"
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi = smooth_voronoi_diagram(voronoi, centerlines, smooth_factor, no_smooth_cl)
        write_polydata(voronoi, voronoi_smoothed_path)

        # Create new surface from the smoothed Voronoi
        surface_smoothed = create_new_surface(voronoi)
        write_polydata(surface_smoothed, surface_smoothed_path)
    elif smooth:
        voronoi = read_polydata(voronoi_smoothed_path)

    return voronoi


def vtk_plane(origin, normal):
    """Returns a vtk box object based on the bounds

    Args:
        origin (list): Center of plane [x, y, z]
        normal (list): Planes normal [x, y, z]

    Returns:
        plane (vtkPlane): A vtkPlane
    """
    plane = vtk.vtkPlane()
    plane.SetOrigin(origin)
    plane.SetNormal(normal)

    return plane


def clip_polydata(surface, cutter=None, value=0):
    """Clip the inpute vtkPolyData object with a cutter function (plane, box, etc)

    Args:
        surface (vtkPolyData): Input vtkPolyData for clipping
        cutter (vtkBox, vtkPlane): Function for cutting the polydata (default None).
        value (float): Distance to the ImplicteFunction or scalar value to clip.

    Returns:
        clipper (vtkPolyData): The clipped surface
    """
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(surface)
    if cutter is None:
        clipper.GenerateClipScalarsOff()
    else:
        clipper.SetClipFunction(cutter)
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(value)
    clipper.Update()

    return clipper.GetOutput(), clipper.GetClippedOutput()


def attach_clipped_regions(surface, clipped, center):
    connectivity = get_connectivity(clipped, mode="All")
    if connectivity.GetNumberOfPoints() == 0:
        return surface
    region_id = get_array("RegionId", connectivity)
    distances = []
    regions = []
    for i in range(int(region_id.max() + 1)):
        regions.append(threshold(connectivity, "RegionId", lower=i - 0.1, upper=i + 0.1, source=0))
        locator = get_locator(regions[-1])
        region_point = regions[-1].GetPoint(locator.FindClosestPoint(center))
        distances.append(distance(region_point, center))

    # Remove the region with the closest distance
    regions.pop(distances.index(min(distances)))

    # Add the other regions back to the surface
    surface = merge_data(regions + [surface])
    surface = clean_surface(surface)
    surface = triangulate_surface(surface)

    return surface


def prepare_surface_output(surface, original_surface, new_centerline, output_filepath,
                           test_merge=False, changed=False, old_centerline=None):
    # Check if the folder for the output exits
    if not path.exists(path.dirname(output_filepath)):
        if path.dirname(output_filepath) != "":
            makedirs(path.dirname(output_filepath))

    # Get planes if outlets of the original surface
    boundary_edges = get_feature_edges(original_surface)
    boundary_connectivity = get_connectivity(boundary_edges)

    vtk_array = boundary_connectivity.GetPointData().GetArray("RegionId")
    vtk_points = boundary_connectivity.GetPoints().GetData()
    region_id = numpy_support.vtk_to_numpy(vtk_array)
    points = numpy_support.vtk_to_numpy(vtk_points)

    centerline = new_centerline if old_centerline is None else old_centerline
    outlets = []
    lines = []
    for i in range(centerline.GetNumberOfLines()):
        lines.append(extract_single_line(centerline, i))
        outlets.append(lines[-1].GetPoint(lines[-1].GetNumberOfPoints() - 1))
    inlet_point = lines[-1].GetPoint(0)

    if changed and old_centerline is None:
        print("WARNING: The changed flag is true, but the old centerline is not provided," + \
              " and the outlet location can therefore not be changed.")

    # Get information from the original geometry
    inlet = False
    for i in range(region_id.max() + 1):
        # Get relevant points
        tmp_points = points[region_id == i]

        # Get normal
        tmp_normal = np.cross(tmp_points[0] - tmp_points[-1],
                              tmp_points[0] - tmp_points[tmp_points.shape[0] // 2])
        normal = tmp_normal / np.sqrt(np.sum(tmp_normal ** 2))

        # Get Center
        center = np.mean(tmp_points, axis=0)

        # Get corresponding centerline to in/outlet
        if np.sqrt(np.sum((np.array(inlet_point) - center) ** 2)) < 0.5:
            line = lines[0]
            line_id = 0
            inlet = True
        else:
            line_id = np.argmin(np.sqrt(np.sum((np.array(outlets) - center) ** 2, axis=1)))
            line = lines[line_id]

        # Set correct direction of normal
        if inlet:
            in_dir = np.array(line.GetPoint(5)) - \
                     np.array(line.GetPoint(0))
        else:
            in_dir = np.array(line.GetPoint(line.GetNumberOfPoints() - 5)) - \
                     np.array(line.GetPoint(line.GetNumberOfPoints() - 1))

        in_dir = in_dir / np.sqrt(np.sum(in_dir ** 2))
        angle = np.arccos(np.dot(in_dir, normal)) * 180 / np.pi
        flipped = True if 90 < angle < 270 else False
        normal = -normal if 90 < angle < 270 else normal
        #print(normal)

        # Mapp the old center and normals to the altered model
        if changed and old_centerline is not None:
            new_line = extract_single_line(new_centerline, line_id)

            # Set correct direction of normal
            if inlet:
                new_outlet = np.array(new_line.GetPoint(0))
                in_dir_new = np.array(new_line.GetPoint(5)) - \
                             new_outlet
                translation = new_outlet - np.array(inlet_point)
            else:
                new_outlet = np.array(new_line.GetPoint(new_line.GetNumberOfPoints() - 1))
                in_dir_new = np.array(new_line.GetPoint(new_line.GetNumberOfPoints() - 5)) - \
                             new_outlet
                translation = new_outlet - np.array(outlets[line_id])

            center += translation
            in_dir_new = in_dir_new / np.sqrt(np.sum(in_dir_new**2))
            in_dir_normal = np.cross(in_dir_new, in_dir)
            dir_angle = np.arccos(np.dot(in_dir, in_dir_new)) * 180 / np.pi

            translation = vtk.vtkTransform()
            translation.RotateWXYZ(-dir_angle, in_dir_normal)
            tmp_normal = normal
            normal = [0, 0, 0]
            translation.TransformNormal(tmp_normal, normal)
            print(" in dir", in_dir)
            print("new dir", in_dir_new)
            print("Old norml", tmp_normal)
            print("New normal", normal)
            print("angle", dir_angle)

        # Set plane
        plane = vtk_plane(center, normal)

        # Clip data (naivly)
        surface, clipped = clip_polydata(surface, plane)

        # Reattach data which should not have been clipped
        surface = attach_clipped_regions(surface, clipped, center)
        inlet = False

    # Perform a 'light' smoothing to obtain a nicer surface
    surface = vmtk_surface_smoother(surface, method="laplace", iterations=100)

    # Clean surface
    surface = clean_surface(surface)
    surface = triangulate_surface(surface)

    # Capped surface
    capped_surface = capp_surface(surface)
    if test_merge:
        check_if_surface_is_merged(capped_surface, new_centerline, output_filepath)

    return surface


def check_if_surface_is_merged(surface, centerlines, output_filepath):
    """
    Clean and check surface for overlapping regions.

    Args:
        surface (vtkPolyData): Surface model.
        centerlines (vtkPolyData): New centerlines.
        output_filepath (str): Filepath of output model.
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
    lines_to_check, _, _ = compute_centerlines(inlet, outlets, None, surface,
                                               resampling=0.1, recompute=True)

    for i in range(centerlines.GetNumberOfLines()):
        line_to_compare = vmtk_centerline_resampling(lines_to_compare[i], length=0.1)
        line_to_check = vmtk_centerline_resampling(extract_single_line(lines_to_check, i), length=0.1)

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
                raise RuntimeError(("\nERROR: Model has most likely overlapping regions." +
                                    " Please check the surface model {} and provide other" +
                                    " parameters for the manipulation or" +
                                    " poly_ball_size.").format(tmp_path))


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

    # Find midpoint and point furthest away
    p1 = P[0]
    p2 = P[1]
    dist = []
    for z in Z:
        d = la.norm(np.cross((z - p1), (z - p2))) / la.norm(p2 - p1)
        dist.append(d)

    d_id = dist.index(max(dist))
    max_dist = max(dist)
    z_m = Z[d_id]

    # Vector from line to Z_max and projection onto plane
    v = (z_m - p2) - (z_m - p2).dot(p2 - p1) * (p2 - p1) / la.norm(p2 - p1) ** 2
    dv = v - v.dot(n) * n

    # Find distances
    dp = p1 + (z_m - p1).dot(p2 - p1) * (p2 - p1) / la.norm(p2 - p1) ** 2
    dpv = dp + dv

    # Move points
    dZ = []
    for i in range(len(Z)):
        dz = np.array(dv) * dist[i] / max_dist * alpha
        dZ.append(Z[i] + dz)

    dx = (dpv - dp) * alpha

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
        # Split points based on orientation to q normal
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
    a, b, c = x[0], x[1], x[1]
    n = np.array([a, b, c])
    n = n / la.norm(n)

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
    a, b, c = dx[0], dx[1], dx[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

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
    a, b, c = dx[0], dx[1], dx[2]
    n = np.array([a, b, c])
    n = n / la.norm(n)

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


def clip_diverging_line(centerline, clip_start_point, clip_end_id):
    """
    Clip the opthamlic artery if present.

    Args:
        centerline (vtkPolyData): Line representing the opthalmic artery centerline.
        clip_start_point (tuple): Point at entrance of opthalmic artery.
        clip_end_id (int): ID of point at end of opthalmic artery.

    Returns:
        patch_eye (vtkPolyData): Voronoi diagram representing opthalmic artery.
    """
    points = [clip_start_point, centerline.GetPoint(clip_end_id)]
    div_points = vtk.vtkPoints()
    for p in points:
        div_points.InsertNextPoint(p)

    patch_cl = CreateParentArteryPatches(centerline, div_points, siphon=True)

    return patch_cl


def get_line_to_change(surface, centerline, region_of_interest, method, region_points,
                       stenosis_length):
    """
    Extract and spline part of centerline
    within the geometry where
    area variations will be performed.

    Args:
        surface (vtkPolyData): Surface model.
        centerline (vtkPolyData): Centerline in geometry.
        region_of_interest (str): Method for setting the region of interest ['manuall' | 'commandline' | 'first_line']
        method (str): Determines which kind of manipulation is performed.
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint.
        stenosis_length (float): Multiplier used to determine the length of the stenosis-affected area.

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
            print("\nPlease select region of interest in the render window.")
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
                elif method == "bend":
                    seed_selector.text = "Press space to select the start and end of the" + \
                                         " bend that you want to manipulate, 'u' to undo.\n"
                seed_selector.Execute()
                stenosis_point_id = seed_selector.GetTargetSeedIds()
                first = False

            region_points = []
            for i in range(stenosis_point_id.GetNumberOfIds()):
                region_points += surface.GetPoint(stenosis_point_id.GetId(i))

        # Get locator
        locator = get_locator(centerline)

        if len(region_points) == 3:
            # Project point onto centerline
            region_points = centerline.GetPoint(locator.FindClosestPoint(region_points))
            point1 = region_points

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
            region_points[:3] = point1
            region_points[3:] = point2

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
                cl_point2 = line.GetPoint(ids2[-1])
                distance1.append(distance(point1, cl_point1))
                distance2.append(distance(point2, cl_point2))

            tol = get_tolerance(centerline) / 10
            total_distance = (np.array(distance1) < tol) * (np.array(distance2) < tol)
            cl_id = np.argmax(total_distance)

            if total_distance[cl_id] == 0:
                raise RuntimeError("The two points provided have to be on the same " + \
                                   " line (from inlet to outlet), and not at two different" + \
                                   " outlets")
            startID = min(ids1[cl_id], ids2[cl_id])
            endID = max(ids1[cl_id], ids2[cl_id])

    # Extract and spline a single line
    line_to_change = extract_single_line(centerline, cl_id, startID=startID, endID=endID)
    line_to_change = spline_centerline(line_to_change, nknots=25, isline=True)

    return line_to_change, region_points


def find_region_of_interest_and_diverging_centerlines(centerlines_complete, region_points):
    """

    Args:
        centerlines_complete (vktPolyData): Complete set of centerlines in geometry.
        region_points (ndarray): Two points determining the region of interest.

    Returns:
        centerlines (vtkPolyData): Centerlines excluding divering lines.
        diverging_centerlines (vtkPolyData): Centerlines diverging from the region of interest.
        region_points (ndarray): Sorted region points.
        region_points_vtk (vtkPoints): Sorted region points as vtkData.
        diverging_ids (list): List of indices where a diverging centerline starts.
    """
    centerlines = []
    diverging_centerlines = []
    p1 = region_points[0]
    p2 = region_points[1]

    # Search for divering centerlines
    tol = get_tolerance(centerlines_complete) * 4
    for i in range(centerlines_complete.GetNumberOfLines()):
        line = extract_single_line(centerlines_complete, i)
        locator = get_locator(line)
        id1 = locator.FindClosestPoint(p1)
        id2 = locator.FindClosestPoint(p2)
        p1_tmp = line.GetPoint(id1)
        p2_tmp = line.GetPoint(id2)
        if distance(p1, p1_tmp) < tol and distance(p2, p2_tmp) < tol:
            centerlines.append(line)
        else:
            diverging_centerlines.append(line)

    # Sort and set clipping points to vtk object
    centerline = centerlines[0]
    locator = get_locator(centerline)
    id1 = locator.FindClosestPoint(region_points[0])
    id2 = locator.FindClosestPoint(region_points[1])
    if id1 > id2:
        region_points = region_points[::-1]
        id1, id2 = id2, id1

    region_points_vtk = vtk.vtkPoints()
    for point in np.asarray(region_points):
        region_points_vtk.InsertNextPoint(point)

    # Find diverging point(s)
    diverging_ids = []
    for line in diverging_centerlines:
        id_end = min([line.GetNumberOfPoints(), centerline.GetNumberOfPoints()])
        for i in np.arange(id1, id_end):
            p_div = np.asarray(line.GetPoint(i))
            p_cl = np.asarray(centerline.GetPoint(i))
            if distance(p_div, p_cl) > tol:
                diverging_ids.append(i)
                break

    centerlines = merge_data(centerlines)
    diverging_centerlines = merge_data(diverging_centerlines) if len(diverging_centerlines) > 0 else None
    return centerlines, diverging_centerlines, region_points, region_points_vtk, diverging_ids


def move_centerlines(patch_cl, dx, p1, p2, diverging_id, diverging_centerlines, direction):
    """

    Args:
        patch_cl (vtkPolyData): Centerlines excluding diverging centerlines.
        dx (ndarray): Direction to move geometry.
        p1 (vtkPolyData): First region point.
        p2: (vtkPolyData): Second region point.
        diverging_id (list): List of index where centerlines diverge from region of interest.
        diverging_centerlines (vtkPolyData): Centerlines which diverge from region of interest.
        direction (str): Manipulation direction parameter.

    Returns:
        centerline (vtkPolyData): Manipulated centerline.
    """
    if diverging_id is not None:
        patch_cl = merge_data([patch_cl, diverging_centerlines])

    numberOfPoints = patch_cl.GetNumberOfPoints()
    numberOfCells = patch_cl.GetNumberOfCells()

    centerline = vtk.vtkPolyData()
    centerlinePoints = vtk.vtkPoints()
    centerlineCellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPoints)

    count = 0
    for i in range(numberOfCells):
        line = extract_single_line(patch_cl, i)
        centerlineCellArray.InsertNextCell(line.GetNumberOfPoints())

        getData = line.GetPointData().GetArray(radiusArrayName).GetTuple1

        locator = get_locator(line)
        id1 = locator.FindClosestPoint(p1)
        if diverging_id is not None and i == (numberOfCells - 1):
            pass
            # Note, reuse id2, idmid from previous loop
        else:
            id2 = locator.FindClosestPoint(p2)
            idmid = int((id1 + id2) * 0.5)

        for p in range(line.GetNumberOfPoints()):
            point = line.GetPoint(p)
            cl_id = locator.FindClosestPoint(point)

            if direction == "horizont":
                if cl_id < id1:
                    dist = dx
                elif id1 <= cl_id < idmid:
                    dist = dx * (idmid ** 2 - cl_id ** 2) / (idmid ** 2 - id1 ** 2)
                elif idmid <= cl_id < (diverging_id - 1) and diverging_id is not None:
                    dist = -dx * (cl_id - idmid) ** 0.5 / (id2 - idmid) ** 0.5
                elif idmid <= cl_id < (id2 - 1):
                    dist = -dx * (cl_id - idmid) ** 0.5 / (id2 - idmid) ** 0.5
                else:
                    if diverging_id is not None and i == (numberOfCells - 1):
                        dist = -dx * (diverging_id - idmid) ** 0.5 / (diverging_id - idmid) ** 0.5
                    else:
                        dist = -dx

            elif direction == "vertical":
                if id1 <= cl_id <= id2:
                    dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2
                elif diverging_id is not None and i == (numberOfCells - 1) and cl_id > id2:
                    cl_id = diverging_id
                    dist = 4 * dx * (cl_id - id1) * (id2 - cl_id) / (id2 - id1) ** 2
                else:
                    dist = 0

            point = np.asarray(point)
            centerlinePoints.InsertNextPoint(point + dist)
            radiusArray.SetTuple1(count, getData(p))
            centerlineCellArray.InsertCellPoint(count)
            count += 1

    centerline.SetPoints(centerlinePoints)
    centerline.SetLines(centerlineCellArray)
    centerline.GetPointData().AddArray(radiusArray)

    return centerline


def split_voronoi_with_centerlines(voronoi, centerline1, centerline2):
    voronoi1 = vtk.vtkPolyData()
    points1 = vtk.vtkPoints()
    cell_array1 = vtk.vtkCellArray()
    radius1 = np.zeros(voronoi.GetNumberOfPoints())
    loc1 = get_locator(centerline1)

    voronoi2 = vtk.vtkPolyData()
    points2 = vtk.vtkPoints()
    cell_array2 = vtk.vtkCellArray()
    radius2 = np.zeros(voronoi.GetNumberOfPoints())
    loc2 = get_locator(centerline2)

    get_radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1

    count1 = 0
    count2 = 0
    for i in range(voronoi.GetNumberOfPoints()):
        point = voronoi.GetPoint(i)
        radius = get_radius(i)
        dist1 = distance(centerline1.GetPoint(loc1.FindClosestPoint(point)), point)
        dist2 = distance(centerline2.GetPoint(loc2.FindClosestPoint(point)), point)

        if dist1 < dist2:
            points1.InsertNextPoint(point)
            radius1[count1] = radius
            cell_array1.InsertNextCell(1)
            cell_array1.InsertCellPoint(count1)
            count1 += 1
        else:
            points2.InsertNextPoint(point)
            radius2[count2] = radius
            cell_array2.InsertNextCell(1)
            cell_array2.InsertCellPoint(count2)
            count2 += 1

    voronoi1.SetPoints(points1)
    voronoi1.SetVerts(cell_array1)
    radius1 = create_vtk_array(radius1[radius1 > 0], radiusArrayName)
    voronoi1.GetPointData().AddArray(radius1)

    voronoi2.SetPoints(points2)
    voronoi2.SetVerts(cell_array2)
    radius2 = create_vtk_array(radius2[radius2 > 0], radiusArrayName)
    voronoi2.GetPointData().AddArray(radius2)

    return voronoi1, voronoi2


def get_clipped_centerline(centerline_relevant_outlets, data):
    line0 = extract_single_line(centerline_relevant_outlets, 0)
    line1 = extract_single_line(centerline_relevant_outlets, 1)
    lines = []
    for i, line in enumerate([line0, line1]):
        loc = get_locator(line)
        tmp_ID_dau = loc.FindClosestPoint(data[i]["end_point"])
        tmp_ID_bif = loc.FindClosestPoint(data["bif"]["end_point"])
        lines.append(extract_single_line(line, 0, startID=tmp_ID_bif, endID=tmp_ID_dau))

    return merge_data(lines)


### The following code is adapted from:
### https://github.com/vmtk/vmtk/tree/master/vmtkApps/CerebralAneurysms/ParentVesselReconstruction
### Written by Marina Piccinelli, and distrubuted within vmtk.
def CreateParentArteryPatches(parentCenterlines, clipPoints, siphon=False, bif=False):
    numberOfDaughterPatches = parentCenterlines.GetNumberOfCells()
    if siphon:
        clipIds, numberOfPatchedCenterlinesPoints = ExtractPatchesIdsSiphon(parentCenterlines, clipPoints)
    else:
        clipIds, numberOfPatchedCenterlinesPoints = ExtractPatchesIds(parentCenterlines, clipPoints)

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

        distance1 = distance(pnt_1, cellLine.GetPoints().GetPoint(ID1))
        distance2 = distance(pnt_2, cellLine.GetPoints().GetPoint(ID2))

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
                                               additionalPoint)

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
                        additionalPoint):
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
        isInside = IsPointInsideInterpolationCylinder(point, cylinderTop, cylinderCenter, cylinderBottom,
                                                      cylinderRadius)

        if isInside == 1:
            maskArray.SetTuple1(i, 1)

    numberOfInterpolationPoints = ComputeNumberOfMaskedPoints(maskArray)

    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfInterpolationPoints)

    count = 0
    for i in range(voronoi.GetNumberOfPoints()):
        value = maskArray.GetTuple1(i)
        if value == 1:
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
            newVoronoiVector = [0.0, 0.0, 0.0]
            newVoronoiPoint = [0.0, 0.0, 0.0]

            transform = vtk.vtkTransform()
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
        arrivalVoronoiPointRadius = voronoiDataset1.GetPointData().GetArray(radiusArrayName).GetTuple1(
            arrivalVoronoiPointId)

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

    if tangentDot < 0.0:
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
                                clippedVoronoi, clippingPoints, bif,
                                cylinder_factor):
    # Extract clipping points
    clippingPointsArray = clippingPoints[1]
    clippingPoints = clippingPoints[0]

    # Copy the voronoi diagram
    completeVoronoiDiagram = vtk.vtkPolyData()
    completeVoronoiDiagram.DeepCopy(clippedVoronoi)

    for j in range(1, 3):
        interpolationCellId = j - 1
        startId = 0
        endId = j

        startCell = vtk.vtkGenericCell()
        patchCenterlines.GetCell(startId, startCell)

        startCellPointId = startCell.GetPointId(startCell.GetNumberOfPoints() - 1)
        startCellPointRadius = patchCenterlines.GetPointData().GetArray(radiusArrayName) \
            .GetTuple1(startCellPointId)
        startCellPointHalfRadius = startCellPointRadius / cylinder_factor

        startInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId,
                                                                                startCellPointId, startCellPointRadius,
                                                                                clippedVoronoi, patchCenterlines)
        startHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(startId,
                                                                                    startCellPointId,
                                                                                    startCellPointHalfRadius,
                                                                                    clippedVoronoi, patchCenterlines)
        endCell = vtk.vtkGenericCell()
        patchCenterlines.GetCell(endId, endCell)

        endCellPointId = endCell.GetPointId(0)
        endCellPointRadius = patchCenterlines.GetPointData().GetArray(radiusArrayName) \
            .GetTuple1(endCellPointId)
        endCellPointHalfRadius = endCellPointRadius / cylinder_factor
        endInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId, endCellPointId,
                                                                              endCellPointRadius, clippedVoronoi,
                                                                              patchCenterlines)
        endHalfInterpolationDataset = ExtractCylindricInterpolationVoronoiDiagram(endId,
                                                                                  endCellPointId,
                                                                                  endCellPointHalfRadius,
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

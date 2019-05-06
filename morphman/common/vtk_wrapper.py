##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

import sys
from os import path

import numpy as np
import vtk

radiusArrayName = 'MaximumInscribedSphereRadius'


def get_number_of_arrays(vtk_polydata):
    """Returns the names and number of arrays for a vtkPolyData object

    Args:
        vtk_polydata (vtkPolyData): object / line to investigate the array.

    Returns:
        count (int): Number of arrays in the line.
        names (list): A list of names of the arrays.
    """
    count = 0
    names = []
    name = 0
    while name is not None:
        name = vtk_polydata.GetPointData().GetArrayName(count)
        if name is not None:
            names.append(name)
            count += 1

    return count, names


def extract_single_line(centerlines, line_id, start_id=0, end_id=None):
    """Extract one line from multiple centerlines.
    If startID and endID is set then only a segment of the centerline is extracted.

    Args:
        centerlines (vtkPolyData): Centerline to extract.
        line_id (int): The line ID to extract.
        start_id (int):
        end_id (int):

    Returns:
        centerline (vtkPolyData): The single line extracted
    """
    cell = vtk.vtkGenericCell()
    centerlines.GetCell(line_id, cell)
    n = cell.GetNumberOfPoints() if end_id is None else end_id + 1

    line = vtk.vtkPolyData()
    cell_array = vtk.vtkCellArray()
    cell_array.InsertNextCell(n - start_id)
    line_points = vtk.vtkPoints()

    arrays = []
    n_, names = get_number_of_arrays(centerlines)
    for i in range(n_):
        tmp = centerlines.GetPointData().GetArray(names[i])
        tmp_comp = tmp.GetNumberOfComponents()
        radius_array = get_vtk_array(names[i], tmp_comp, n - start_id)
        arrays.append(radius_array)

    point_array = []
    for i in range(n_):
        point_array.append(centerlines.GetPointData().GetArray(names[i]))

    count = 0
    for i in range(start_id, n):
        cell_array.InsertCellPoint(count)
        line_points.InsertNextPoint(cell.GetPoints().GetPoint(i))

        for j in range(n_):
            num = point_array[j].GetNumberOfComponents()
            if num == 1:
                tmp = point_array[j].GetTuple1(i)
                arrays[j].SetTuple1(count, tmp)
            elif num == 2:
                tmp = point_array[j].GetTuple2(i)
                arrays[j].SetTuple2(count, tmp[0], tmp[1])
            elif num == 3:
                tmp = point_array[j].GetTuple3(i)
                arrays[j].SetTuple3(count, tmp[0], tmp[1], tmp[2])
            elif num == 9:
                tmp = point_array[j].GetTuple9(i)
                arrays[j].SetTuple9(count, tmp[0], tmp[1], tmp[2], tmp[3], tmp[4],
                                    tmp[5], tmp[6], tmp[7], tmp[8])
        count += 1

    line.SetPoints(line_points)
    line.SetLines(cell_array)
    for j in range(n_):
        line.GetPointData().AddArray(arrays[j])

    return line


def read_polydata(filename, datatype=None):
    """
    Load the given file, and return a vtkPolyData object for it.

    Args:
        filename (str): Path to input file.
        datatype (str): Additional parameter for vtkIdList objects.

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
        # Read header
        with open(filename) as myfile:
            head = "".join(list(islice(myfile, 10))).lower()

        # Set reader based on header
        if "unstructured_grid" in head:
            reader = vtk.vtkUnstructuredGridReader()
        elif "structured_grid" in head:
            reader = vtk.vtkStructuredGridReader()
        elif "rectilinear_grid" in head:
            reader = vtk.vtkRectilinearGridReader()
        elif "structured_points" in head:
            reader = vtk.vtkStructuredPointsReader()
        elif "polydata" in head:
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
    elif fileType == "np" and datatype == "vtkIdList":
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
    polydata = reader.GetOutput()

    return polydata


def write_polydata(input_data, filename, datatype=None, file_type="ascii"):
    """
    Write the given input data based on the file name extension.

    Args:
        file_type (string): Filetype of output
        input_data (vtkSTL/vtkPolyData/vtkXMLStructured/
                    vtkXMLRectilinear/vtkXMLPolydata/vtkXMLUnstructured/
                    vtkXMLImage/Tecplot): Input data.
        filename (str): Save path location.
        datatype (str): Additional parameter for vtkIdList objects.
    """
    # Check filename format
    fileType = filename.split(".")[-1]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')

    # Get writer
    if fileType == 'stl':
        writer = vtk.vtkSTLWriter()

    elif fileType == 'vtk':
        # Set reader based on data type
        if isinstance(input_data, vtk.vtkUnstructuredGrid):
            writer = vtk.vtkUnstructuredGridWriter()
        elif isinstance(input_data, vtk.vtkStructuredGrid):
            writer = vtk.vtkStructuredGridWriter()
        elif isinstance(input_data, vtk.vtkRectilinearGrid):
            writer = vtk.vtkRectilinearGridWriter()
        elif isinstance(input_data, vtk.vtkStructuredPoints) or \
                isinstance(input_data, vtk.vtkImageData):
            writer = vtk.vtkStructuredPointsWriter()
        elif isinstance(input_data, vtk.vtkPolyData):
            writer = vtk.vtkPolyDataWriter()

        if file_type.lower() == "ascii":
            writer.SetFileType(1)
        elif file_type.lower() == "binary":
            writer.SetFileType(0)
        else:
            raise ValueError("Invalid file type, can only be ascii or binary")

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
    elif fileType == "np" and datatype == "vtkIdList":
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


def vtk_merge_polydata(inputs):
    """
    Appends one or more polygonal
    datates together into a single
    polygonal dataset.

    Args:
        inputs (list): List of vtkPolyData objects.

    Returns:
        merged_data (vtkPolyData): Single polygonal dataset.
    """
    append_filter = vtk.vtkAppendPolyData()
    for input_ in inputs:
        append_filter.AddInputData(input_)
    append_filter.Update()
    merged_data = append_filter.GetOutput()

    return merged_data


def get_point_data_array(array_name, line, k=1):
    """
    Get data array from centerline object (Point data).

    Args:
        array_name (str): Name of array.
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
        data_array = line.GetPointData().GetArray(array_name).GetTuple1
    elif k == 2:
        data_array = line.GetPointData().GetArray(array_name).GetTuple2
    elif k == 3:
        data_array = line.GetPointData().GetArray(array_name).GetTuple3

    for i in range(line.GetNumberOfPoints()):
        array[i, :] = data_array(i)

    return array


def write_vtk_points(points, filename):
    """
    Writes input points to file.

    Args:
        points (vtkPolyData): Point data.
        filename (str): Save location.
    """
    point_set = vtk.vtkPolyData()
    cell_array = vtk.vtkCellArray()

    for i in range(points.GetNumberOfPoints()):
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i)

    point_set.SetPoints(points)
    point_set.SetVerts(cell_array)

    write_polydata(point_set, filename)


def vtk_clean_polydata(surface):
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
    surface_cleaner = vtk.vtkCleanPolyData()
    surface_cleaner.SetInputData(surface)
    surface_cleaner.Update()
    cleaned_surface = surface_cleaner.GetOutput()

    return cleaned_surface


def vtk_compute_connectivity(surface, mode="All", closest_point=None, show_color_regions=True,
                             mark_visited_points=False):
    """Wrapper of vtkPolyDataConnectivityFilter. Compute connectivity.

    Args:
        show_color_regions (bool): Turn on/off the coloring of connected regions.
        mark_visited_points (bool): Specify whether to record input point ids that appear in the output.
        surface (vtkPolyData): Input surface data.
        mode (str): Type of connectivity filter.
        closest_point (list): Point to be used for mode='Closest'
    """
    connectivity = vtk.vtkPolyDataConnectivityFilter()
    connectivity.SetInputData(surface)

    # Mark each region with "RegionId"
    if mode == "All":
        connectivity.SetExtractionModeToAllRegions()
    elif mode == "Largest":
        connectivity.SetExtractionModeToLargestRegion()
    elif mode == "Closest":
        if closest_point is None:
            print("ERROR: point not set for extracting closest region")
            sys.exit(0)
        connectivity.SetExtractionModeToClosestPointRegion()
        connectivity.SetClosestPoint(closest_point)

    if show_color_regions:
        connectivity.ColorRegionsOn()

    if mark_visited_points:
        connectivity.MarkVisitedPointIdsOn()

    connectivity.Update()
    output = connectivity.GetOutput()

    return output


def vtk_convert_unstructured_grid_to_polydata(unstructured_grid):
    """Wrapper for vtkGeometryFilter, which converts an unstructured grid into a polydata.

    Args:
        unstructured_grid (vtkUnstructuredGrid): An unstructured grid.

    Returns:
        surface (vtkPolyData): A vtkPolyData object from the unstrutured grid.
    """
    # Convert unstructured grid to polydata
    geo_filter = vtk.vtkGeometryFilter()
    geo_filter.SetInputData(unstructured_grid)
    geo_filter.Update()
    polydata = geo_filter.GetOutput()

    return polydata


def vtk_compute_threshold(surface, name, lower=0, upper=1, threshold_type="between", source=1):
    """Wrapper for vtkThreshold. Extract a section of a surface given a criteria.

    Args:
        surface (vtkPolyData): The input data to be extracted.
        name (str): Name of scalar array.
        lower (float): Lower bound.
        upper (float): Upper bound.
        threshold_type (str): Type of threshold (lower, upper, between)
        source (int): PointData or CellData.

    Returns:
        surface (vtkPolyData): The extracted surface based on the lower and upper limit.
    """
    # source = 1 uses cell data as input
    # source = 0 uses point data as input

    # Apply threshold
    vtk_threshold = vtk.vtkThreshold()
    vtk_threshold.SetInputData(surface)
    if threshold_type == "between":
        vtk_threshold.ThresholdBetween(lower, upper)
    elif threshold_type == "lower":
        vtk_threshold.ThresholdByLower(lower)
    elif threshold_type == "upper":
        vtk_threshold.ThresholdByUpper(upper)
    else:
        print((("%s is not a threshold type. Pleace chose from: upper, lower" +
                ", or between") % threshold_type))
        sys.exit(0)

    vtk_threshold.SetInputArrayToProcess(0, 0, 0, source, name)
    vtk_threshold.Update()
    surface = vtk_threshold.GetOutput()

    # Convert to polydata
    surface = vtk_convert_unstructured_grid_to_polydata(surface)

    return surface


def vtk_extract_feature_edges(polydata, compute_feature_edges=False, compute_boundary_edges=True,
                              compute_non_manifold_edges=False):
    """Wrapper for vtkFeatureedges. Extracts the edges of the cells that are open.

    Args:
        compute_non_manifold_edges (bool): Turn on/off the extraction of non-manifold edges.
        compute_boundary_edges (bool): Turn on/off the extraction of boundary edges.
        compute_feature_edges (bool): Turn on/off the extraction of feature edges.
        polydata (vtkPolyData): surface to extract the openings from.

    Returns:
        feature_edges (vtkPolyData): The boundary edges of the surface.
    """
    feature_edges = vtk.vtkFeatureEdges()
    if compute_feature_edges:
        feature_edges.FeatureEdgesOn()
    else:
        feature_edges.FeatureEdgesOff()
    if compute_boundary_edges:
        feature_edges.BoundaryEdgesOn()
    else:
        feature_edges.BoundaryEdgesOff()
    if compute_non_manifold_edges:
        feature_edges.NonManifoldEdgesOn()
    else:
        feature_edges.NonManifoldEdgesOff()
    feature_edges.SetInputData(polydata)
    feature_edges.Update()

    return feature_edges.GetOutput()


def get_vtk_array(name, comp, num):
    """An empty vtkDoubleArray.

    Args:
        name (str): Name of array.
        comp (int): Number of components
        num (int): Number of tuples.

    Returns:
        array (vtkDoubleArray): An empty vtk array.
    """
    array = vtk.vtkDoubleArray()
    array.SetNumberOfComponents(comp)
    array.SetNumberOfTuples(num)
    for i in range(comp):
        array.FillComponent(i, 0.0)
    array.SetName(name)

    return array


def get_vtk_cell_locator(surface):
    """Wrapper for vtkCellLocator

    Args:
        surface (vtkPolyData): input surface

    Returns:
        return (vtkCellLocator): Cell locator of the input surface.
    """
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(surface)
    locator.BuildLocator()

    return locator


def create_vtk_array(values, name, k=1):
    """Given a set of numpy values, and a name of the array create vtk array

    Args:
        values (numpy.ndarray): List of the values.
        name (str): Name of the array.
        k (int): Length of tuple.

    Returns:
        vtk_array (vtkPointArray): vtk point array
    """
    vtk_array = get_vtk_array(name, k, values.shape[0])

    if k == 1:
        for i in range(values.shape[0]):
            vtk_array.SetTuple1(i, values[i])
    elif k == 2:
        for i in range(values.shape[0]):
            vtk_array.SetTuple2(i, values[i, 0], values[i, 1])
    elif k == 3:
        for i in range(values.shape[0]):
            vtk_array.SetTuple3(i, values[i, 0], values[i, 1], values[i, 2])
    elif k == 9:
        for i in range(values.shape[0]):
            vtk_array.SetTuple9(i, values[i, 0], values[i, 1], values[i, 2],
                                values[i, 3], values[i, 4], values[i, 5],
                                values[i, 6], values[i, 7], values[i, 8])

    return vtk_array


def move_past_sphere(centerline, center, r, start, step=-1, stop=0, scale_factor=0.8):
    """Moves a point along the centerline until it as outside the a sphere with radius (r)
    and a center (center).

    Args:
        centerline (vtkPolyData): Centerline to move along.
        center (list): point list of the center of the sphere
        r (float): the radius of a sphere
        start (int): id of the point along the centerline where to start.
        step (int): direction along the centerline.
        stop (int): ID along centerline, for when to stop searching.
        scale_factor (float): Scale the radius with this factor.

    Returns:
        tmp_point (list): The first point on the centerline outside the sphere
        r (float): minimal inscribed sphere radius at the new point.
        i (int): the centerline ID at the new point.
    """
    # Create the minimal inscribed sphere
    misr_sphere = vtk.vtkSphere()
    misr_sphere.SetCenter(center)
    misr_sphere.SetRadius(r * scale_factor)
    tmp_point = [0.0, 0.0, 0.0]

    # Go the length of one MISR backwards
    for i in range(start, stop, step):
        value = misr_sphere.EvaluateFunction(centerline.GetPoint(i))
        if value >= 0.0:
            tmp_point = centerline.GetPoint(i)
            break

    r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

    return tmp_point, r, i


def vtk_point_locator(centerline):
    """Wrapper for vtkStaticPointLocator.

    Args:
        centerline (vtkPolyData): Input vtkPolyData.

    Returns:
        locator (vtkStaticPointLocator): Point locator of the input surface.
    """
    locator = vtk.vtkStaticPointLocator()
    locator.SetDataSet(centerline)
    locator.BuildLocator()

    return locator


def vtk_triangulate_surface(surface, pass_lines=False, pass_verts=False):
    """Wrapper for vtkTriangleFilter.

    Args:
        pass_lines (bool): Turn on/off passing lines through filter. Default is On
        pass_verts (bool): Turn on/off passing vertices through filter. Default is On
        surface (vtkPolyData): Surface to triangulate.

    Returns:
        surface (vtkPolyData): Triangulated surface.
    """
    surface_triangulator = vtk.vtkTriangleFilter()
    surface_triangulator.SetInputData(surface)
    if not pass_lines:
        surface_triangulator.PassLinesOff()
    if not pass_verts:
        surface_triangulator.PassVertsOff()
    surface_triangulator.Update()

    return surface_triangulator.GetOutput()


def vtk_compute_mass_properties(surface, compute_surface_area=True, compute_volume=False):
    """
    Calculate the volume from the given polydata

    Args:
        compute_volume (bool): Compute surface volume if True
        compute_surface_area (bool): Compute surface area if True
        surface (vtkPolyData): Surface to compute are off

    Returns:
        area (float): Area of the input surface
    Returns:
        volume (float): Volume of the input surface
    """
    mass = vtk.vtkMassProperties()
    mass.SetInputData(surface)

    if compute_surface_area:
        return mass.GetSurfaceArea()

    if compute_volume:
        return mass.GetVolume()


def vtk_marching_cube(modeller, compute_normals=False, compute_scalars=False, compute_gradients=False):
    """Wrapper for vtkMarchingCube

    Args:
        modeller (vtkPolyballModeller): Modeller of a surface model
        compute_normals (bool): Set/Get the computation of normals.
        compute_scalars (bool): Set/Get the computation of scalars.
        compute_gradients (bool): Set/Get the computation of gradients.

    Returns:
        vtkMarchingCube: Isosurface generated from surface
    """
    marching_cube = vtk.vtkMarchingCubes()
    if compute_normals:
        marching_cube.ComputeNormalsOn()
    if compute_scalars:
        marching_cube.ComputeScalarsOn()
    if compute_gradients:
        marching_cube.ComputeGradientsOn()

    marching_cube.SetInputData(modeller.GetOutput())
    marching_cube.SetValue(0, 0.0)
    marching_cube.Update()

    return marching_cube


def vtk_compute_normal_gradients(cell_normals, use_faster_approximation=False):
    """
    Compute gradients of the normals

    Args:
        cell_normals (vtkPolyData): Surface to compute normals on
        use_faster_approximation (bool): Use a less accurate algorithm that performs fewer calculations, but faster.
    """
    gradient_filter = vtk.vtkGradientFilter()
    gradient_filter.SetInputData(cell_normals)
    gradient_filter.SetInputArrayToProcess(0, 0, 0, 1, "Normals")
    if use_faster_approximation:
        gradient_filter.FasterApproximationOn()

    gradient_filter.Update()
    gradients = gradient_filter.GetOutput()

    return gradients


def vtk_compute_polydata_normals(surface, compute_point_normals=False, compute_cell_normals=False):
    """ Wrapper for vtkPolyDataNormals

    Args:
        surface (vtkPolyData): Surface model
        compute_point_normals (bool): Turn on/off the computation of point normals.
        compute_cell_normals (bool): Turn on/off the computation of cell normals.

    Returns:
        vtkPolyData: Cell normals of surface model
    """
    normal_generator = vtk.vtkPolyDataNormals()
    normal_generator.SetInputData(surface)
    if compute_point_normals:
        normal_generator.ComputePointNormalsOn()
    else:
        normal_generator.ComputePointNormalsOff()
    if compute_cell_normals:
        normal_generator.ComputeCellNormalsOn()
    else:
        normal_generator.ComputeCellNormalsOff()

    normal_generator.Update()
    cell_normals = normal_generator.GetOutput()

    return cell_normals


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


def vtk_clip_polydata(surface, cutter=None, value=0, get_inside_out=False, generate_clip_scalars=False):
    """Clip the inpute vtkPolyData object with a cutter function (plane, box, etc)

    Args:
        generate_clip_scalars (bool): If True, output scalar values will be interpolated from implicit function values.
        get_inside_out (bool): Get inside out, default is False
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
    if get_inside_out:
        clipper.InsideOutOn()
    if generate_clip_scalars and cutter is not None:
        clipper.GenerateClipScalarsOn()
    clipper.GenerateClippedOutputOn()
    clipper.SetValue(value)
    clipper.Update()

    return clipper.GetOutput(), clipper.GetClippedOutput()

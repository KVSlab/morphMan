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


def write_polydata(input_data, filename, datatype=None):
    """
    Write the given input data based on the file name extension.

    Args:
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


def vtk_append_polydata(inputs):
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

def write_points(points, filename):
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



def vtk_compute_connectivity(surface, mode="All", closest_point=None):
    """Wrapper of vtkPolyDataConnectivityFilter. Compute connectivity.

    Args:
        surface (vtkPolyData): Input surface data.
        mode (str): Type of connectivity filter.
        closest_point (list): Point to be used for mode='Closest'"""
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
    connectivity.ColorRegionsOn()
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


def vtk_extract_feature_edges(polydata):
    """Wrapper for vtkFeatureedges. Extracts the edges of the cells that are open.

    Args:
        polydata (vtkPolyData): surface to extract the openings from.

    Returns:
        feature_edges (vtkPolyData): The boundary edges of the surface.
    """
    feature_edges = vtk.vtkFeatureEdges()
    feature_edges.FeatureEdgesOff()
    feature_edges.BoundaryEdgesOn()
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


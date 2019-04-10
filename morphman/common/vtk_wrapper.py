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
    append_filter = vtk.vtkAppendPolyData()
    for input_ in inputs:
        append_filter.AddInputData(input_)
    append_filter.Update()
    merged_data = append_filter.GetOutput()

    return merged_data


def get_array(array_name, line, k=1):
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


def get_array_cell(array_name, line, k=1):
    """
    Get data array from centerline object (Cell data).

    Args:
        array_name (str): Name of array.
        line (vtkPolyData): Centerline object.
        k (int): Dimension.

    Returns:
        array (ndarray): Array containing data points.
    """
    # vtk_array = line.GetCellData().GetArray(arrayName)
    # array = numpy_support.vtk_to_numpy(vtk_array)
    array = np.zeros((line.GetNumberOfCells(), k))
    if k == 1:
        data_array = line.GetCellData().GetArray(array_name).GetTuple1
    elif k == 2:
        data_array = line.GetCellData().GetArray(array_name).GetTuple2
    elif k == 3:
        data_array = line.GetCellData().GetArray(array_name).GetTuple3
    elif k == 9:
        data_array = line.GetCellData().GetArray(array_name).GetTuple9

    for i in range(line.GetNumberOfCells()):
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


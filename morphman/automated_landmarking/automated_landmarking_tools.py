##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from scipy.ndimage.filters import gaussian_filter
from scipy.signal import argrelextrema
from vmtk import vmtkrenderer
from vtk.numpy_interface import dataset_adapter as dsa

# Local import
from morphman.common import *


def get_maximum_coronal_coordinate(coronal_coordinate, curvilinear_coordinates):
    """
    Find the maximum coronal coordinate and return its ID along the centerline

    Args:
        coronal_coordinate (ndarray): Coronal coordinates along the centerline
        curvilinear_coordinates (ndarray): Curvilinear coordinates along the centerline

    Returns:
        max_coronal_coordinate_id (int): ID representing the maximum coronal coordinate along the centerline
    """
    value_index = coronal_coordinate[argrelextrema(coronal_coordinate, np.less_equal)[0]].min()
    max_coronal_coordinate_id = np.array(coronal_coordinate.tolist().index(value_index))
    search_tolerance = 30

    if abs(curvilinear_coordinates[max_coronal_coordinate_id] - curvilinear_coordinates[-1]) > search_tolerance:
        print("-- Sanity check failed, checking for maximums")

        value_index = coronal_coordinate[argrelextrema(coronal_coordinate, np.greater_equal)[0]].max()
        max_coronal_coordinate_id = np.array(coronal_coordinate.tolist().index(value_index))
        if abs(curvilinear_coordinates[max_coronal_coordinate_id] - curvilinear_coordinates[-1]) > search_tolerance:
            print("-- Sanity check failed, no anterior bend in model. Exiting.")
            sys.exit(1)

    return max_coronal_coordinate_id


def get_centerline_coordinates(line, curvilinear_coordinates):
    """
    Get coordinates along the centerline from the curvilinear coordinates

    Args:
        line (vtkPolyData): Centerline
        curvilinear_coordinates (ndarray): Curvilinear coordinates

    Returns:
        (dict): Dictionary containing centerline coordinates
    """
    x, y, z = dsa.WrapDataObject(line).GetPoints()[:curvilinear_coordinates.shape[0]].T

    return dict(x=x, y=y, z=z)


def orient_centerline(ica_centerline):
    """
    Check the orientation of the axial coordinate and
    reverse centerline if path of centerline is in distal direction.

    Args:
        ica_centerline (vtkPolyData): Centerline to check

    Returns:
        ica_centerline (vtkPolyData): Centerline, possibly reversed
    """
    curvilinear_coordinates = get_curvilinear_coordinate(ica_centerline)
    coordinates = get_centerline_coordinates(ica_centerline, curvilinear_coordinates)

    # Iterate through coordinates and find axial coordinate
    for _, coordinate in coordinates.items():
        coordinate = gaussian_filter(coordinate, 25)
        id0 = np.where(coordinate == np.min(coordinate))[0][0]
        id1 = np.where(coordinate == np.max(coordinate))[0][0]
        if id0 == 0 and id1 == (len(coordinate) - 1):
            return ica_centerline

        if id1 == 0 and id0 == (len(coordinate) - 1):
            return reverse_centerline(ica_centerline)

    return ica_centerline


def spline_centerline_and_compute_geometric_features(line, smooth, nknots):
    """
    Compute attributes and geometric parameters of input
    centerline, using B-splines (SciPy).

    Args:
        line (vtkPolyData): Centerline data.
        smooth (bool): Smooth centerline with VMTK if True.
        nknots (int): Number of knots for B-splines.

    Returns:
        line (vtkPolyData): Splined centerline.
    Returns:
        max_point_ids (ndarray): Array of max curvature values
    Returns:
        min_point_ids (ndarray): Array of min curvature values
    """
    data = dsa.WrapDataObject(line).GetPoints()
    curv_coor = get_curvilinear_coordinate(line)

    t = np.linspace(curv_coor[0], curv_coor[-1], nknots + 2)[1:-1]

    fx = splrep(curv_coor, data[:, 0], k=4, t=t)
    fy = splrep(curv_coor, data[:, 1], k=4, t=t)
    fz = splrep(curv_coor, data[:, 2], k=4, t=t)

    data[:, 0] = splev(curv_coor, fx)
    data[:, 1] = splev(curv_coor, fy)
    data[:, 2] = splev(curv_coor, fz)

    header = ["X", "Y", "Z"]
    line = convert_numpy_data_to_polydata(data, header)

    # Let vmtk compute curve attributes
    line = vmtk_compute_centerline_attributes(line)
    line = vmtk_compute_geometric_features(line, smooth)

    # Compute curvature from the 'exact' spline to get a robust way of
    # finding max / min points on the centerline
    dlsfx = splev(curv_coor, fx, der=1)
    dlsfy = splev(curv_coor, fy, der=1)
    dlsfz = splev(curv_coor, fz, der=1)

    ddlsfx = splev(curv_coor, fx, der=2)
    ddlsfy = splev(curv_coor, fy, der=2)
    ddlsfz = splev(curv_coor, fz, der=2)

    c1xc2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
    c1xc2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
    c1xc2_3 = ddlsfy * dlsfx - ddlsfx * dlsfy

    curvature_ = np.sqrt(c1xc2_1 ** 2 + c1xc2_2 ** 2 + c1xc2_3 ** 2) / \
                 (dlsfx ** 2 + dlsfy ** 2 + dlsfz ** 2) ** 1.5

    max_point_ids = list(argrelextrema(curvature_, np.greater)[0])
    min_point_ids = list(argrelextrema(curvature_, np.less)[0])

    locator = get_vtk_point_locator(line)

    min_points = [[data[:, 0][i], data[:, 1][i], data[:, 2][i]] for i in min_point_ids]
    max_points = [[data[:, 0][i], data[:, 1][i], data[:, 2][i]] for i in max_point_ids]
    min_point_ids = []
    max_point_ids = []

    for point_min, point_max in zip(min_points, max_points):
        min_point_ids.append(locator.FindClosestPoint(point_min))
        max_point_ids.append(locator.FindClosestPoint(point_max))

    curvature = get_point_data_array("Curvature", line)
    line = get_k1k2_basis(curvature, line)

    length = get_curvilinear_coordinate(line)
    dddlsfx = splev(length, fx, der=3)
    dddlsfy = splev(length, fy, der=3)
    dddlsfz = splev(length, fz, der=3)

    torsion_spline = (dddlsfx * c1xc2_1 + dddlsfy * c1xc2_2 + dddlsfz * c1xc2_3) / \
                     (c1xc2_1 ** 2 + c1xc2_2 ** 2 + c1xc2_3 ** 2)
    torsion_array = create_vtk_array(torsion_spline, "Torsion")
    line.GetPointData().AddArray(torsion_array)

    curvature_ = np.sqrt(c1xc2_1 ** 2 + c1xc2_2 ** 2 + c1xc2_3 ** 2) / \
                 (dlsfx ** 2 + dlsfy ** 2 + dlsfz ** 2) ** 1.5
    curvature_[0] = curvature[0]
    curvature_[-1] = curvature[-1]

    curvature_array = create_vtk_array(curvature_, "Curvature")
    line.GetPointData().AddArray(curvature_array)

    return line, max_point_ids, min_point_ids


def map_landmarks(landmarks, centerline, algorithm):
    """
    Takes new landmarks and original centerline,
    mapping each landmark interface to the original centerline.
    Filters away duplicate landmarks

    Args:
        landmarks (dict): Contains landmarks.
        centerline (vtkPolyData): Original centerline.
        algorithm (str): Landmarking algorith, bogunovic or piccinelli

    Returns:
        mapped_landmarks (dict): Contains landmarks mapped to centerline.
    """
    mapped_landmarks = {}
    landmark_ids = []
    locator = get_vtk_point_locator(centerline)
    k = 1
    for key in landmarks:
        landmark = landmarks[key]
        landmark_id = locator.FindClosestPoint(landmark)

        if algorithm in ["piccinelli", "kjeldsberg"]:
            key = "C%s" if algorithm == "kjeldsberg" else "bend%s"
            if landmark_id in landmark_ids:
                continue  # Skip for duplicate landmarking point
            else:
                landmark_ids.append(landmark_id)
                landmark_mapped = centerline.GetPoint(landmark_id)
                mapped_landmarks[key % k] = landmark_mapped
                k += 1

        if algorithm == "bogunovic":
            landmark_mapped = centerline.GetPoint(landmark_id)
            landmarks[key] = landmark_mapped
            mapped_landmarks = landmarks

    return mapped_landmarks


def create_particles(base_path, algorithm, method, step=1):
    """
    Create a file with points where bends are located and
    remove points from manifest

    Args:
        base_path (str): Case location.
        algorithm (str): Name of landmarking algorithm.
        method (str): Method used for computing curvature.
    """

    info_filepath = base_path + "_info.json"
    filename_all_landmarks = base_path + "_landmark_%s_%s_step%s.particles" % (algorithm, method, step)
    print("Saving all landmarks to: %s" % filename_all_landmarks)

    output_all = open(filename_all_landmarks, "w")
    with open(info_filepath, ) as landmarks_json:
        landmarked_points = json.load(landmarks_json)

    for key in landmarked_points:
        if algorithm == "bogunovic":
            if key in ["anterior_posterior", "posterior_inferior", "inferior_end", "superior_anterior"]:
                p = landmarked_points[key]
                point = "%s %s %s" % (p[0], p[1], p[2])
                output_all.write(point + "\n")

        elif algorithm == "piccinelli":
            if key[:4] == "bend":
                p = landmarked_points[key]
                point = "%s %s %s" % (p[0], p[1], p[2])
                output_all.write(point + "\n")

        elif algorithm == "kjeldsberg":
            if key[0] == "C":
                p = landmarked_points[key]
                point = "%s %s %s" % (p[0], p[1], p[2])
                output_all.write(point + "\n")

    output_all.close()


def mark_diverging_arteries(centerline_complete, ica_centerline):
    print("\nPlease select the diverging arteries in the render window.")

    first = True
    while True:
        if not first:
            print("Please provide only one or no points, try again")

        # Select point on surface
        seed_selector = vmtkPickPointSeedSelector()
        seed_selector.NoPointSelection = True
        seed_selector.SetSurface(centerline_complete)
        seed_selector.text = "Press space to select a point along the centerline, corresponding " + \
                             "to the Ophthalmic Artery OR press 'q' to skip. Press 'u' to undo.\n"

        seed_selector.Execute()
        ophthalmic_id = seed_selector.GetTargetSeedIds()
        first = False
        if ophthalmic_id.GetNumberOfIds() in [0, 1]:
            break

    first = True
    while True:
        if not first:
            print("Please provide only one or no points, try again")

        # Select point on surface
        seed_selector = vmtkPickPointSeedSelector()
        seed_selector.NoPointSelection = True
        seed_selector.SetSurface(centerline_complete)
        seed_selector.text = "Press space to select a point along the centerline, corresponding " + \
                             "to the Posterior Communicating Artery OR press 'q' to skip. Press 'u' to undo.\n"

        seed_selector.Execute()
        p_com_id = seed_selector.GetTargetSeedIds()
        first = False
        if p_com_id.GetNumberOfIds() in [0, 1]:
            break

    tol = 4 * get_centerline_tolerance(ica_centerline)

    if ophthalmic_id.GetNumberOfIds() == 1:
        ophthalmic = centerline_complete.GetPoint(ophthalmic_id.GetId(0))
        ophthalmic_line = get_closest_line(ophthalmic, centerline_complete)
        ophthalmic_div_id = get_diverging_point_id(ophthalmic_line, ica_centerline, tol)
    else:
        ophthalmic_div_id = None

    if p_com_id.GetNumberOfIds() == 1:
        p_com = centerline_complete.GetPoint(p_com_id.GetId(0))
        p_com_line = get_closest_line(p_com, centerline_complete)
        p_com_div_id = get_diverging_point_id(p_com_line, ica_centerline, tol)
    else:
        p_com_div_id = None

    classify_with_diverging_arteries = True if ophthalmic_div_id or p_com_div_id is not None else False

    return classify_with_diverging_arteries, ophthalmic_div_id, p_com_div_id


def get_closest_line(reference_point, centerline_complete):
    closest_line = None
    dist = 1E9
    for i in range(centerline_complete.GetNumberOfLines()):
        tmp_line = extract_single_line(centerline_complete, i)
        locator = get_vtk_point_locator(tmp_line)
        point_id = locator.FindClosestPoint(reference_point)
        point = tmp_line.GetPoint(point_id)

        tmp_dist = la.norm(np.asarray(point) - np.asarray(reference_point))
        if tmp_dist < dist:
            dist = tmp_dist
            closest_line = tmp_line

    return closest_line


def check_for_diverging_centerlines(centerline_complete, line):
    classify_using_diverging_centerlines = False
    region_points = [np.asarray(line.GetPoint(0)), np.asarray(line.GetPoint(line.GetNumberOfPoints() - 5))]
    centerline_diverging = get_region_of_interest_and_diverging_centerlines(centerline_complete, region_points)

    p_com_a_id = None
    ophthalmic_id = None

    if len(centerline_diverging[-1]) >= 2:
        diverging_ids = sorted(centerline_diverging[-1])
        ophthamlic = line.GetPoint(diverging_ids[0])
        p_com_a = line.GetPoint(diverging_ids[1])
        ica_locator = get_vtk_point_locator(line)
        ophthalmic_id = ica_locator.FindClosestPoint(ophthamlic)
        p_com_a_id = ica_locator.FindClosestPoint(p_com_a)
        classify_using_diverging_centerlines = True

    return classify_using_diverging_centerlines, ophthalmic_id, p_com_a_id


def visualize_landmarks(landmarks, centerline, algorithm, surface_path):
    """
    Visualize the interfaces between the resulting landmarks (points)

    Args:
        landmarks (ndarray): Points defining landmarks / interfaces between segments
        centerline (vtkPolyData): Centerline which has been landmarked
        algorithm (str): Name of landmarking algorithm
        surface_path (str): Location of surface model
    """

    vtkPoints = vtk.vtkPoints()
    landmarkSpheres = []
    for key, point in landmarks.items():
        sphere = vtk.vtkSphereSource()
        sphere.SetRadius(1)
        sphere.SetCenter(point[0], point[1], point[2])
        sphere.Update()
        landmarkSpheres.append(sphere.GetOutput())
        vtkPoints.InsertNextPoint(point)

    landmarkData = vtk_merge_polydata(landmarkSpheres)

    surface = None
    if path.exists(surface_path):
        surface = read_polydata(surface_path)

    vtkNameArray = vtk.vtkStringArray()
    vtkNameArray.SetNumberOfValues(7)

    if algorithm == "piccinelli":
        start = len(landmarks) - 1
        stop = -1
        interfaces = ["Bend %s-%s" % (i, i + 1) for i in range(len(landmarks))]
    elif algorithm == "bogunovic":
        start = 3
        stop = -1
        interfaces = [
            "Anterior-Posterior",
            "Posterior-Inferior",
            "Inferior (Start)",
            "Superior-Anterior",
        ]
    elif algorithm == "kjeldsberg":
        start = 6
        stop = 7 - len(landmarks) - 1
        interfaces = ["C%s-C%s" % (i - 1, i) for i in range(1, len(landmarks) + 1)]

    for i in range(start, stop, -1):
        if i == 0:
            if algorithm == "kjeldsberg":
                vtkNameArray.SetValue(i, 'C1 (Start)')
            elif algorithm == "piccinelli":
                vtkNameArray.SetValue(i, 'Bend 1 (Start)')
            elif algorithm == "bogunovic":
                vtkNameArray.SetValue(i, interfaces[0])
        else:
            vtkNameArray.SetValue(i, interfaces[i])

    # Create the dataset with landmarking points
    labeledPolydata = vtk.vtkPolyData()
    labeledPolydata.SetPoints(vtkPoints)
    labeledPolydata.GetPointData().AddArray(vtkNameArray)

    centerlineMapper = vtk.vtkPolyDataMapper()
    centerlineMapper.SetInputData(centerline)
    centerlineActor = vtk.vtkActor()
    centerlineActor.SetMapper(centerlineMapper)
    centerlineActor.GetProperty().SetLineWidth(10)
    centerlineActor.GetProperty().SetColor(1, 0, 0)

    landmarkMapper = vtk.vtkPolyDataMapper()
    landmarkMapper.SetInputData(landmarkData)
    landmarkActor = vtk.vtkActor()
    landmarkActor.SetMapper(landmarkMapper)

    labelsMapper = vtk.vtkLabeledDataMapper()
    labelsMapper.SetInputData(labeledPolydata)
    labelsMapper.SetLabelModeToLabelFieldData()
    labelsMapper.GetLabelTextProperty().SetColor(1, 1, 1)
    labelsMapper.GetLabelTextProperty().SetFontSize(35)
    labelsMapper.GetLabelTextProperty().ItalicOff()

    labelsActor = vtk.vtkActor2D()
    labelsActor.SetMapper(labelsMapper)

    vmtkRenderer = vmtkrenderer.vmtkRenderer()
    vmtkRenderer.Initialize()
    vmtkRenderer.Renderer.AddActor(centerlineActor)
    vmtkRenderer.Renderer.AddActor(labelsActor)
    vmtkRenderer.Renderer.AddActor(landmarkActor)

    if surface is not None:
        surfaceMapper = vtk.vtkPolyDataMapper()
        surfaceMapper.SetInputData(surface)
        surfaceActor = vtk.vtkActor()
        surfaceActor.SetMapper(surfaceMapper)
        surfaceActor.GetProperty().SetOpacity(0.3)
        vmtkRenderer.Renderer.AddActor(surfaceActor)

    vmtkRenderer.Render()

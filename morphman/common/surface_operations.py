def create_new_surface(complete_voronoi_diagram, poly_ball_size=[120, 120, 120]):
    """
    Envelops an input voronoi diagram
    into a new surface model at a
    given resolution determined by
    the poly_ball_size.

    Args:
        complete_voronoi_diagram (vtkPolyData): Voronoi diagram
        poly_ball_size (list): List of dimensional resolution of output model

    Returns:
        envelope (vtkPolyData): Enveloped surface model.
    """
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    modeller.SetInputData(complete_voronoi_diagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(poly_ball_size)
    modeller.Update()

    # Write the new surface
    marching_cube = vtk.vtkMarchingCubes()
    marching_cube.SetInputData(modeller.GetOutput())
    marching_cube.SetValue(0, 0.0)
    marching_cube.Update()
    envelope = marching_cube.GetOutput()

    return envelope

def compute_circleness(surface):
    """Compute the area ratio betwen minimum circle and the maximum circle.

    Args:
        surface (vtkPolyData): Boundary edges of an opening

    Returns:
        circleness (float): Area ratio
        center (list): Center of the opening.
    """
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

def check_if_surface_is_capped(surface):
    """Checks if the surface is closed, and how many openings there are.

    Args:
        surface (vtkPolyData): Surface to be checked

    Returns:
        open (boolean): Open or closed surface
        number (int): Number of integer
    """
    # Get boundary cells
    cells = get_feature_edges(surface)
    if cells.GetNumberOfCells() == 0:
        return True, 0
    else:
        outlets = get_connectivity(cells, mode="All")
        number = get_array("RegionId", outlets).max()
        return number == 0, int(number)

def uncapp_surface(surface, gradients_limit=0.15, area_limit=0.3, circleness_limit=3):
    """
    A rule-based method for removing endcapps on a surface. The method considers the
    gradient of the normals, the size of the region, and how similar it is to a circle.

    Args:
        surface (vtkPolyData): Surface to be uncapped.
        gradients_limit (float): Upper limit for gradients of normals.
        area_limit (float): Lower limit of the area.
        circleness_limit (float): Upper limit of the circleness.

    Returns:
        surface (vtkPolyData): The uncapped surface.

    """

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
                          threshold_type="between", source=1)

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
                                 upper=(i + limit), threshold_type="between", source=0))
        circ, center = compute_circleness(regions[-1])
        circleness.append(circ)
        centers_edge.append(center)
        area.append(compute_area(regions[-1]))

    # Only keep outlets with circleness < circleness_limit and area > area_limit
    circleness_ids = np.where(np.array(circleness) < circleness_limit)
    region_ids = np.where(np.array(area) > area_limit)
    regions = [regions[i] for i in region_ids[0] if i in circleness_ids[0]]
    centers_edge = [centers_edge[i] for i in region_ids[0] if i in circleness_ids[0]]

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
            cell_id = vtk.mutable(0)
            sub_id = vtk.mutable(0)
            dist = vtk.mutable(0)
            locator.FindClosestPoint(tmp_center, p, cell_id, sub_id, dist)
            mark_outlets.SetTuple(cell_id, [1])

    surface.GetCellData().AddArray(mark_outlets)

    # Remove the outlets from the original surface
    uncapped_surface = threshold(surface, "outlets", lower=0, upper=0.5,
                                 threshold_type="between", source=1)

    # Check if some cells where not marked
    remove = True
    while remove:
        locator = get_locator_cell(uncapped_surface)
        mark_outlets = create_vtk_array(np.zeros(uncapped_surface.GetNumberOfCells()), "outlets", k=1)
        remove = False
        for center in centers_edge:
            locator.FindClosestPoint(center, p, cell_id, sub_id, dist)
            if dist < 0.01:
                remove = True
                mark_outlets.SetTuple(cell_id, [1])

        uncapped_surface.GetCellData().AddArray(mark_outlets)

        if remove:
            uncapped_surface = threshold(uncapped_surface, "outlets", lower=0,
                                         upper=0.5, threshold_type="between", source=1)

    return uncapped_surface

def clipp_capped_surface(surface, centerlines, clipspheres=0):
    """A method for clipping a capped outlets. The branches will be clipped some distance
    from the outlets.

    Args:
        surface (vtkPolyData): Surface to clipp
        centerlines (vtkPolyData): Centerlines to mark the in and outlets.
        clipspheres (float): Number of end point spheres

    Returns:
        surface (vtkPolyData): Clipped surface
    """
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


def compute_area(surface):
    """
    Compute area of polydata

    Args:
        surface (vtkPolyData): Surface to compute are off

    Returns:
        area (float): Area of the input surface
    """
    mass = vtk.vtkMassProperties()
    mass.SetInputData(surface)


    return mass.GetSurfaceArea()

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
    print("-- Please select the two relevant outlets in the interactive window.")
    seed_selector = vmtkPickPointSeedSelector()
    seed_selector.SetSurface(triangulated_surface)
    seed_selector.text = "Please select the two relevant outlets, \'u\' to undo\n"
    seed_selector.Execute()

    point_seed_ids = seed_selector.GetTargetSeedIds()
    get_point = surface.GetPoints().GetPoint
    points = [list(get_point(point_seed_ids.GetId(i))) for i in range(point_seed_ids.GetNumberOfIds())]
    info = {}

    if dir_path is not None:
        for i in range(len(points)):
            info["relevant_outlet_%d" % i] = points[i]
        write_parameters(info, dir_path)

    return points


def triangulate_surface(surface):
    """Wrapper for vtkTriangleFilter.

    Args:
        surface (vtkPolyData): Surface to triangulate.

    Returns:
        surface (vtkPolyData): Triangulated surface.
    """
    surface_triangulator = vtk.vtkTriangleFilter()
    surface_triangulator.SetInputData(surface)
    surface_triangulator.PassLinesOff()
    surface_triangulator.PassVertsOff()
    surface_triangulator.Update()

    return surface_triangulator.GetOutput()

def get_centers(surface, base_path, flowext=False):
    """Get the centers of the inlet and outlets.

    Args:
        surface (vtkPolyData): An open surface.
        base_path (str): Path to the case file.
        flowext (bool): Turn on/off flow extension.

    Returns:
        inlet (list): A flatt list with the point of the inlet
        outlet (list): A flatt list with the points of all the outlets.
    """
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


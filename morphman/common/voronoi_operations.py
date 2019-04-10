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

    new_voronoi = voronoi.VoronoiDiagram
    return new_voronoi

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
    number_of_points = voronoi.GetNumberOfPoints()
    thresholds = get_array(radiusArrayName, centerlines) * (1 - smoothing_factor)

    # Do not smooth inlet and outlets, set threshold to -1
    start = 0
    end = 0
    for i in range(centerlines.GetNumberOfLines()):
        line = extract_single_line(centerlines, i)
        length = get_curvilinear_coordinate(line)
        end_ = line.GetNumberOfPoints() - 1
        end += end_

        # Point buffer start
        end_id = end_ - np.argmin(np.abs(-(length - length.max()) - thresholds[end]))
        start_id = np.argmin(np.abs(length - thresholds[start]))

        thresholds[start:start + start_id] = -1
        thresholds[end - end_id:end] = -1
        start += end_ + 1
        end += 1

    locator = get_locator(centerlines)
    if no_smooth_cl is not None:
        no_locator = get_locator(no_smooth_cl)

    smoothed_diagram = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cell_array = vtk.vtkCellArray()
    radius_array_numpy = np.zeros(number_of_points)

    count = 0
    for i in range(number_of_points):
        point = voronoi.GetPoint(i)
        radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
        id_ = locator.FindClosestPoint(point)
        cl_point = centerlines.GetPoint(id_)

        if distance(point, cl_point) > 2 * thresholds[id_] / (1 - smoothing_factor):
            points.InsertNextPoint(point)
            cell_array.InsertNextCell(1)
            cell_array.InsertCellPoint(count)
            radius_array_numpy[count] = radius
            count += 1

        elif no_smooth_cl is not None:
            dist1 = distance(point, centerlines.GetPoint(id_))
            id_1 = no_locator.FindClosestPoint(point)
            dist2 = distance(point, no_smooth_cl.GetPoint(id_1))

            if dist2 < dist1:
                points.InsertNextPoint(point)
                cell_array.InsertNextCell(1)
                cell_array.InsertCellPoint(count)
                radius_array_numpy[count] = radius
                count += 1
            else:
                if radius >= thresholds[id_]:
                    points.InsertNextPoint(point)
                    cell_array.InsertNextCell(1)
                    cell_array.InsertCellPoint(count)
                    radius_array_numpy[count] = radius
                    count += 1
        else:
            if radius >= thresholds[id_]:
                points.InsertNextPoint(point)
                cell_array.InsertNextCell(1)
                cell_array.InsertCellPoint(count)
                radius_array_numpy[count] = radius
                count += 1

        radius_array = get_vtk_array(radiusArrayName, 1, count)

    for i in range(count):
        radius_array.SetTuple1(i, radius_array_numpy[i])

    smoothed_diagram.SetPoints(points)
    smoothed_diagram.SetVerts(cell_array)
    smoothed_diagram.GetPointData().AddArray(radius_array)

    return smoothed_diagram




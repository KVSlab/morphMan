##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from morphman.common.centerline_operations import *


def remove_distant_voronoi_points(voronoi, centerline):
    """Take a voronoi diagram and a centerline remove points that are far away.

    Args:
        voronoi (vtkPolyData): Voronoi data.
        centerline (vtkPolyData): centerline.

    Returns:
        voronoi (vtkPolyData): Voronoi diagram without the extreme points
    """
    n = voronoi.GetNumberOfPoints()
    new_voronoi = vtk.vtkPolyData()
    cell_array = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    radius = np.zeros(n)

    locator = vtk_point_locator(centerline)
    radius_array_data = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1
    limit = radius_array_data(0)
    limit = limit * 10

    count = 0
    for i in range(n):
        point = voronoi.GetPoint(i)
        cl_point_id = locator.FindClosestPoint(point)
        cl_point = centerline.GetPoint(cl_point_id)
        dist = get_distance(point, cl_point)
        if dist / 3 > radius_array_data(i) or radius_array_data(i) > limit:
            count += 1
            continue

        points.InsertNextPoint(point)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i - count)
        value = radius_array_data(i)
        radius[i - count] = value

    print("Removed %s points from the voronoi diagram" % count)

    radius_array = get_vtk_array(radiusArrayName, 1, n - count)
    for i in range(n - count):
        radius_array.SetTuple(i, [float(radius[i])])

    new_voronoi.SetPoints(points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)

    return new_voronoi


def smooth_voronoi_diagram(voronoi, centerlines, smoothing_factor, no_smooth_cl=None):
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
    thresholds = get_point_data_array(radiusArrayName, centerlines) * (1 - smoothing_factor)

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

    locator = vtk_point_locator(centerlines)
    if no_smooth_cl is not None:
        no_locator = vtk_point_locator(no_smooth_cl)

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

        if get_distance(point, cl_point) > 2 * thresholds[id_] / (1 - smoothing_factor):
            points.InsertNextPoint(point)
            cell_array.InsertNextCell(1)
            cell_array.InsertCellPoint(count)
            radius_array_numpy[count] = radius
            count += 1

        elif no_smooth_cl is not None:
            dist1 = get_distance(point, centerlines.GetPoint(id_))
            id_1 = no_locator.FindClosestPoint(point)
            dist2 = get_distance(point, no_smooth_cl.GetPoint(id_1))

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
    modeller = vmtk_polyball_modeller(complete_voronoi_diagram, poly_ball_size)

    # Write the new surface
    marching_cube = vtk_marching_cube(modeller)
    envelope = marching_cube.GetOutput()

    return envelope


def get_split_voronoi_diagram(voronoi, centerlines):
    """Given two centerlines, and a Voronoi diagram, return two Voronoi diagrams based on
    the distance of the two centerlines.

    Args:
        voronoi (vtkPolyData): Input Voronoi diagram
        centerlines (list): A list of centerlines (vtkPolyData). An entery could
                            alternativly be None as well, the corresponding voronoi
                            diagram would then be None as well.

    Returns
        voronoi2 (list): A list of Voronoi diagrams closest to each centerline.
    """
    n = len(centerlines)
    centerline1 = [centerlines[i] for i in range(n) if centerlines[i] is not None]
    voronoi1 = [vtk.vtkPolyData() for i in range(n) if centerlines[i] is not None]
    points1 = [vtk.vtkPoints() for i in range(n) if centerlines[i] is not None]
    cell_array1 = [vtk.vtkCellArray() for i in range(n) if centerlines[i] is not None]
    radius1 = [np.zeros(voronoi.GetNumberOfPoints()) for i in range(n) if centerlines[i] is not None]
    loc1 = [vtk_point_locator(centerlines[i]) for i in range(n) if centerlines[i] is not None]
    count1 = [0 for i in range(n) if centerlines[i] is not None]

    n1 = len(centerline1)

    get_radius = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1

    for i in range(voronoi.GetNumberOfPoints()):
        dists = []
        point = voronoi.GetPoint(i)
        radius = get_radius(i)
        for i in range(n1):
            dists.append(get_distance(centerline1[i].GetPoint(loc1[i].FindClosestPoint(point)), point))

        index = dists.index(min(dists))

        points1[index].InsertNextPoint(point)
        radius1[index][count1[index]] = radius
        cell_array1[index].InsertNextCell(1)
        cell_array1[index].InsertCellPoint(count1[index])
        count1[index] += 1

    for i in range(n1):
        voronoi1[i].SetPoints(points1[i])
        voronoi1[i].SetVerts(cell_array1[i])
        tmp_radius1 = create_vtk_array(radius1[i][radius1[i] > 0], radiusArrayName)
        voronoi1[i].GetPointData().AddArray(tmp_radius1)

    if n1 != n:
        voronoi2 = []
        for i in range(n):
            if centerlines[i] is None:
                voronoi2.append(None)
            else:
                voronoi2.append(voronoi1[i])
    else:
        voronoi2 = voronoi1

    return voronoi2

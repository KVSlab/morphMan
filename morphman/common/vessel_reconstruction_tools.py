import math

from scipy.interpolate import splrep, splev

from morphman.common.common import get_distance
from morphman.common.vmtk_wrapper import *
from morphman.common.vtk_wrapper import *


### The following code is adapted from:
### https://github.com/vmtk/vmtk/tree/master/vmtkApps/CerebralAneurysms/ParentVesselReconstruction
### Written by Marina Piccinelli, and distrubuted within vmtk.
def create_parent_artery_patches(parentCenterlines, clipPoints, siphon=False, bif=False):
    """Clip out a segment of the centerline, and create new centerlines with new end and
    starting points.

    Args:
        parentCenterlines (vtkPolyData): Original centerline
        clipPoints (vtkPoints): The points where to clip the centerline.
        siphon (bool): On/off clipping a siphon
        bif (bool): On/off bifurcation.

    Returns:
        centerline (vtkPolyData): New centerline without the segment.
    """
    numberOfDaughterPatches = parentCenterlines.GetNumberOfCells()
    if siphon:
        clipIds, numberOfPatchedCenterlinesPoints = extract_patches_ids_siphon(parentCenterlines,
                                                                               clipPoints)
    else:
        clipIds, numberOfPatchedCenterlinesPoints = extract_patches_ids(parentCenterlines,
                                                                        clipPoints)

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


def extract_patches_ids_siphon(parentCl, clipPts, clipped=False):
    """For each clipping points (clipPts) extract the cooresponding ID for each line in
    the centerline. (This is for the siphon, see extract_patches_ids as well.)

    Args:
        parentCl (vtkPolyData):
        clipPts (vtkPoints):
        clipped (bool):

    Returns:
        clipIds (list): A list of IDs.
        numberOfPoints (int): Total number of points.
    """
    clipIds = []
    numberOfPoints = 0

    upstreamPoint = clipPts.GetPoint(0)
    downstreamPoint = clipPts.GetPoint(1)

    for j in range(parentCl.GetNumberOfCells()):
        cellLine = extract_single_line(parentCl, j)
        locator = vtk_point_locator(cellLine)

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


def extract_patches_ids(parentCl, clipPts):
    """For each clipping points (clipPts) extract the cooresponding ID for each line in
    the centerline.

    Args:
        parentCl (vtkPolyData):
        clipPts (vtkPoints):

    Returns:
        clipIds (list): A list of IDs.
        numberOfPoints (int): Total number of points.
    """
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
        locator = vtk_point_locator(cellLine)

        if j == 0 and N == 3:
            upstreamId = locator.FindClosestPoint(commonPoint)
            clipIds.append(upstreamId)
            numberOfPoints += upstreamId + 1

        ID1 = locator.FindClosestPoint(pnt_1)
        ID2 = locator.FindClosestPoint(pnt_2)

        distance1 = get_distance(pnt_1, cellLine.GetPoints().GetPoint(ID1))
        distance2 = get_distance(pnt_2, cellLine.GetPoints().GetPoint(ID2))

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


def interpolate_patch_centerlines(patchCenterlines, parentCenterlines,
                                  additionalPoint, lower, version, tension=0,
                                  continuity=0):
    """Interpolate new centerlines between end and starting points. Given
    additionalPoiint, lower, and version, then number and method for interpolation varies.

    Args:
        patchCenterlines (vtkPolyData): Clipped centerline.
        parentCenterlines (vtkPolyData): The original centerline.
        additionalPoint (vtkPoints): Additional point to interpolate through.
        lower (str): None / 'lower' / 'bif' to indicate how to interpolate.
        version (bool): Method for interpolation.
        tension (float): Variable for the Kochanek spline
        continuity (float): Variable for the Kochanek spline

    Returns:
        centerline (vtkPolyData): The new centerline, including the new interpolated
        segment.
    """
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
            splinePoints = interpolate_spline(startingCell, endingCell, additionalPoint)
        else:
            splinePoints = interpolate_two_cells(startingCell, endingCell,
                                                 numberOfInterpolationPoints,
                                                 additionalPointIds[i], additionalPoint,
                                                 tension, continuity)

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


def interpolate_spline(startCell, endCell, additionalPoint):
    """Interpolate between two lines using splrep from scipy, potentially with an
    additional point (additionalPoint).

    Args:
        startCell (vtkPolyData): Start line
        endCell (tkPolyData): End line
        additionalPoint (list): A list with the coordinates to the additional point.

    Returns:
        centerline (vtkPolyData): The new interpolated centerline.
    """
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
        curv_coor[i + 1] = curv_coor[i] + get_distance(points[i], points[i + 1])

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


def interpolate_two_cells(startCell, endCell, numberOfSplinePoints, additionalPointId,
                          additionalPoint, tension, continuitiy):
    """Interpolate between two lines using vtkCardinalSpline from vtk, potentially with an
    additional point (additionalPoint).

    Args:
        startCell (vtkPolyData): Start line
        endCell (tkPolyData): End line
        numberOfSplinePoints (int): Number of spline point.
        additionalPointId (int): Id of the additional point.
        additionalPoint (list): A list with the coordinates to the additional point.
        tension (float): Variable for the Kochanek spline
        continuity (float): Variable for the Kochanek spline

    Returns:
        centerline (vtkPolyData): The new interpolated centerline.
    """
    points = vtk.vtkPoints()
    # xspline = vtk.vtkCardinalSpline()
    # yspline = vtk.vtkCardinalSpline()
    # zspline = vtk.vtkCardinalSpline()
    xspline = vtk.vtkKochanekSpline()
    yspline = vtk.vtkKochanekSpline()
    zspline = vtk.vtkKochanekSpline()

    for s in [xspline, yspline, zspline]:
        s.SetDefaultTension(tension)
        s.SetDefaultContinuity(continuitiy)

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


def extract_cylindric_interpolation_voronoi_diagram(cellId, pointId, cylinderRadius,
                                                    voronoi, centerlines,
                                                    interpolationHalfSize=3):
    """Extract the voronoi diagram within a cylinder to be used for extrapolation.

    Args:
        cellId (int): LineId of the centerline.
        pointId (int): Point Id of where to extract the cylinder.
        cylinderRadius (float): The radius of the cylinder.
        voronoi (vtkPolyData): The voronoi diagram to extract cylinder from.
        centerlines (vtkPolyData): Centerline corresponding to the Voronoi diagram.

    Returns:
        interpolationDataset (vtkPolyData): The extracted cylinder from the Voronoi
        diagram.
    """
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
        isInside = is_point_inside_interpolation_cylinder(point, cylinderTop,
                                                          cylinderCenter, cylinderBottom,
                                                          cylinderRadius)

        if isInside == 1:
            maskArray.SetTuple1(i, 1)

    numberOfInterpolationPoints = compute_number_of_masked_points(maskArray)

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


def is_point_inside_interpolation_cylinder(x, t, c, b, r):
    """Check if a (Voronoi) point is inside a cylinder.

    Args:
        x (list): Point to check.
        t (list): Top of the cylinder.
        c (list): Center of the cylinder.
        b (list): Bottom of the cylinder.
        r (float): Radius of the cylinder.

    Returns:
        inside (bool): True if inside, False if outside.
    """
    halfheigth = get_distance(b, t) / 2

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


def compute_number_of_masked_points(data_array):
    """Count number of '1' in the data array.

    Args:
        data_array (vtkIntArray): Array to check.

    Returns:
        number_of_points (int): Number of '1' in array.
    """
    number_of_points = 0
    for i in range(data_array.GetNumberOfTuples()):
        value = data_array.GetTuple1(i)
        if value == 1:
            number_of_points += 1

    return number_of_points


def voronoi_diagram_interpolation(interpolationcellid, id0, id1, voronoiDataset0,
                                  voronoiDataset1, centerlines, step,
                                  clippingPoints):
    """Given two Voronoi datasets interpolate the data sets along the centerline.

    Args:
        interpolationcellid (int): LineID of the centerline
        id0 (int): Start ID.
        id1 (int): Stop ID.
        voronoiDataset0 (vtkPolyData): First Voronoi dataset.
        voronoiDataset1 (vtkPolyData): Second Voronoi dataset.
        centerlines (vtkPolyData): Centerline to interpolate along.
        step (int): Direction to interpolate
        clippingPoints (vtkPoints): Location of clipping points.

    Returns:
        finalNewVoronoiPoints (vtkPoints): New points to the Voronoi diagram.
        finalRadiusArray (vtkDoubelArray): Array to hold the radius for each point.
    """
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

        centerlinePointLocator = vtk_point_locator(cellLine)

        closestPointId = centerlinePointLocator.FindClosestPoint(voronoiPoint)
        closestPoint = cellLine.GetPoint(closestPointId)

        voronoiVector = [0.0, 0.0, 0.0]
        voronoiVector[0] = voronoiPoint[0] - closestPoint[0]
        voronoiVector[1] = voronoiPoint[1] - closestPoint[1]
        voronoiVector[2] = voronoiPoint[2] - closestPoint[2]
        voronoiVectorNorm = vtk.vtkMath.Norm(voronoiVector)
        rotationAngle = compute_voronoi_vector_to_centerline_angle(closestPointId, voronoiVector, cellLine)

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

        voronoiPointLocator = vtk_point_locator(voronoiDataset1)

        arrivalVoronoiPointId = voronoiPointLocator.FindClosestPoint(lastPTPoint)
        arrivalVoronoiPoint = voronoiDataset1.GetPoint(arrivalVoronoiPointId)
        arrivalVoronoiPointRadius = voronoiDataset1.GetPointData().GetArray(radiusArrayName).GetTuple1(
            arrivalVoronoiPointId)

        arrivalCenterlinePointLocator = vtk_point_locator(cellLine)

        arrivalCenterlineClosestPointId = arrivalCenterlinePointLocator.FindClosestPoint(arrivalVoronoiPoint)
        arrivalCenterlineClosestPoint = cellLine.GetPoint(arrivalCenterlineClosestPointId)

        arrivalVoronoiVector = [0.0, 0.0, 0.0]
        arrivalVoronoiVector[0] = arrivalVoronoiPoint[0] - arrivalCenterlineClosestPoint[0]
        arrivalVoronoiVector[1] = arrivalVoronoiPoint[1] - arrivalCenterlineClosestPoint[1]
        arrivalVoronoiVector[2] = arrivalVoronoiPoint[2] - arrivalCenterlineClosestPoint[2]
        arrivalVoronoiVectorNorm = vtk.vtkMath.Norm(arrivalVoronoiVector)
        radiusArray = compute_spline(voronoiPointRadius, arrivalVoronoiPointRadius, numberOfPTPoints)
        vectorNormArray = compute_spline(voronoiVectorNorm, arrivalVoronoiVectorNorm, numberOfPTPoints)

        pointsToGap = (gapStartId - closestPointId) * step
        if pointsToGap < 0 or PTPoints.GetNumberOfPoints() <= pointsToGap:
            continue

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


def compute_voronoi_vector_to_centerline_angle(pointId, vector, centerline):
    """Compute angle between the normal compontent of the centerline, and the parallel
    transport normal.

    Args:
        pointId (int): Id along centerline of interest.
        vector (list): Vector
        centerline (vtkPolyData): Centerline

    Returns:
        alpha (float): Angle
    """

    point0 = centerline.GetPoint(pointId)
    point1 = centerline.GetPoint(pointId + 1)
    point2 = centerline.GetPoint(pointId - 1)

    tangent = [0.0, 0.0, 0.0]
    for i in range(3):
        tangent[i] += point1[i] - point0[i]
        tangent[i] += point0[i] - point2[i]

    ptnnormal = centerline.GetPointData().GetArray(parallelTransportNormalsArrayName).GetTuple3(pointId)
    alpha = compute_angle_between_vectors(ptnnormal, tangent, vector)

    return alpha


def normalize(vector):
    """Normalize a vector to unit length.

    Args:
        vector (list/tuple/numpy.ndarray): Array to be normalized

    Return:
        vector (numpy.ndarray): Normalized vector.
    """
    length = np.sqrt(np.sum(np.asarray(vector) ** 2))
    if length == 0:
        return np.asarray(vector)
    else:
        return np.asarray(vector) / length


def compute_angle_between_vectors(normal, tangent, vector):
    """Compute the angle between the vector and the normal component.

    Args:
        normal (list): Normal component to the centerline.
        tangent (list): Tangent to the centerline.
        vector (list): Vector to compute the angle of.

    Returns:
        angle (float): Angle in degrees.
    """
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


def compute_spline(startValue, endValue, numberOfPoints):
    """Create a vtkDoubleArray which is a spline between startValue, mean(startValue,
    endValue), and endValue.

    Args:
        startValue (float): The start value.
        endValue (float): The end value.
        numberOfPoints (int): The number of points.
    """
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


def insert_new_voronoi_points(oldDataset, newPoints, newArray):
    """Insert new points into an unconnected vtkPolyData object.

    Args:
        oldDataset (vtkPolyData): object to insert new points.
        newPoints (vtkPoints): New points to insert into the oldDataset.
        newArray (vtkDoubleArray): Array corresponding to the points.

    Returns:
        newDataset (vtkPolyData): the oldDataset, but with the new points.
    """
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

    count = 0
    for i in range(numberOfNewPoints):
        point = [0.0, 0.0, 0.0]
        newPoints.GetPoint(i, point)

        if np.sum(np.isinf(point) + np.isnan(point)) > 0:
            continue

        if np.sum(np.abs(np.array(point)) > 10000) > 0:
            continue

        points.InsertNextPoint(point)
        cellArray.InsertNextCell(1)
        cellArray.InsertCellPoint(numberOfDatasetPoints + count)

        value = newArray.GetTuple1(i)
        radiusArray.SetTuple1(numberOfDatasetPoints + count, value)
        count += 1

    newDataset.SetPoints(points)
    newDataset.SetVerts(cellArray)
    newDataset.GetPointData().AddArray(radiusArray)

    return newDataset


def interpolate_voronoi_diagram(interpolatedCenterlines, patchCenterlines,
                                clippedVoronoi, clippingPoints, bif,
                                cylinder_factor):
    """Extrapolate/interpolate a the Voronoi diagram along new centerlines. This is the
    core of the algorithm to reconstruct a bifurcation in manipulate_branches.py.

    Args:
        interpolatedCenterlines (vtkPolyData): Centerlines which has been interpolated.
        patchCenterlines (vtkPolyData): Centerlines without the interpolated patch.
        clippedVoronoi (vtkPolyData): Clipped Voronoi diagram.
        clippingPoints (vtkPoints): Points at where the centerline and Voronoi diagram
        where clipped.
        bif (list): List of extra centerlines to extrapolate along.
        cylinder_factor (float): Size of cylinder to extract as the basis of the
        extrapolation of the Voronoi diagram.

    Returns:
        completeVoronoiDiagram (vtkPolyData): The modified Voronoi diagram.
    """
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

        startInterpolationDataset = extract_cylindric_interpolation_voronoi_diagram(startId,
                                                                                    startCellPointId,
                                                                                    startCellPointRadius,
                                                                                    clippedVoronoi, patchCenterlines)
        startHalfInterpolationDataset = extract_cylindric_interpolation_voronoi_diagram(startId,
                                                                                        startCellPointId,
                                                                                        startCellPointHalfRadius,
                                                                                        clippedVoronoi,
                                                                                        patchCenterlines)
        endCell = vtk.vtkGenericCell()
        patchCenterlines.GetCell(endId, endCell)

        endCellPointId = endCell.GetPointId(0)
        endCellPointRadius = patchCenterlines.GetPointData().GetArray(radiusArrayName) \
            .GetTuple1(endCellPointId)
        endCellPointHalfRadius = endCellPointRadius / cylinder_factor
        endInterpolationDataset = extract_cylindric_interpolation_voronoi_diagram(endId, endCellPointId,
                                                                                  endCellPointRadius, clippedVoronoi,
                                                                                  patchCenterlines)
        endHalfInterpolationDataset = extract_cylindric_interpolation_voronoi_diagram(endId,
                                                                                      endCellPointId,
                                                                                      endCellPointHalfRadius,
                                                                                      clippedVoronoi, patchCenterlines)

        # Find and insert new points
        newVoronoiPoints, newVoronoiPointsMISR = voronoi_diagram_interpolation(interpolationCellId,
                                                                               startId, endId,
                                                                               startInterpolationDataset,
                                                                               endHalfInterpolationDataset,
                                                                               interpolatedCenterlines, 1,
                                                                               clippingPoints)
        completeVoronoiDiagram = insert_new_voronoi_points(completeVoronoiDiagram, newVoronoiPoints,
                                                           newVoronoiPointsMISR)

        completeVoronoiDiagram.GetNumberOfPoints()
        newVoronoiPoints, newVoronoiPointsMISR = voronoi_diagram_interpolation(interpolationCellId,
                                                                               endId, startId, endInterpolationDataset,
                                                                               startHalfInterpolationDataset,
                                                                               interpolatedCenterlines, -1,
                                                                               clippingPoints)
        completeVoronoiDiagram = insert_new_voronoi_points(completeVoronoiDiagram, newVoronoiPoints,
                                                           newVoronoiPointsMISR)

    if bif is not []:
        for i in range(len(bif) - 1):
            interpolationCellId = 0
            startId_ = 0
            endId_ = 1
            bif_clipped = bif[-1]
            bif_ = bif[i]

            startCell = extract_single_line(bif_clipped, startId_)
            locator = vtk_point_locator(startCell)
            startPoint = clippingPoints.GetPoint(1)
            startId = locator.FindClosestPoint(startPoint)
            startR = startCell.GetPointData().GetArray(radiusArrayName).GetTuple1(startId)
            startRHalf = startR / cylinder_factor

            endCell = extract_single_line(bif_clipped, endId_)
            locator = vtk_point_locator(endCell)
            endPoint = endCell.GetPoint(0)
            endId = locator.FindClosestPoint(endPoint)
            endR = endCell.GetPointData().GetArray(radiusArrayName).GetTuple1(endId)
            endRHalf = endR / cylinder_factor

            id1, id2 = get_start_ids(clippingPointsArray, bif_clipped)

            startInterpolationDataset = extract_cylindric_interpolation_voronoi_diagram(startId_,
                                                                                        startId, startR,
                                                                                        clippedVoronoi, startCell)
            startHalfInterpolationDataset = extract_cylindric_interpolation_voronoi_diagram(startId_,
                                                                                            startId, startRHalf,
                                                                                            clippedVoronoi, startCell)

            endInterpolationDataset = extract_cylindric_interpolation_voronoi_diagram(endId_, endId,
                                                                                      endR, clippedVoronoi,
                                                                                      endCell)
            endHalfInterpolationDataset = extract_cylindric_interpolation_voronoi_diagram(endId_, endId,
                                                                                          endRHalf, clippedVoronoi,
                                                                                          endCell)
            newVoronoiPoints, newVoronoiPointsMISR = voronoi_diagram_interpolation(interpolationCellId,
                                                                                   id1, id2, startInterpolationDataset,
                                                                                   endHalfInterpolationDataset,
                                                                                   bif_, 1, clippingPoints)

            completeVoronoiDiagram = insert_new_voronoi_points(completeVoronoiDiagram, newVoronoiPoints,
                                                               newVoronoiPointsMISR)
            newVoronoiPoints, newVoronoiPointsMISR = voronoi_diagram_interpolation(interpolationCellId,
                                                                                   id2, id1, endInterpolationDataset,
                                                                                   startHalfInterpolationDataset,
                                                                                   bif_, -1, clippingPoints)

            completeVoronoiDiagram = insert_new_voronoi_points(completeVoronoiDiagram, newVoronoiPoints,
                                                               newVoronoiPointsMISR)

    return completeVoronoiDiagram


def get_start_ids(points, line):
    """Sort the points according the the distance from the line.

    Args:
        points (list): Nested-list with points.
        line (vtkPolyData): Centerline of interest.

    Return:
        order1 (int): Proximal point
        order2 (int): Distal point
    """
    p1 = points[1]
    p2 = points[2]
    locator = vtk_point_locator(line)
    id1 = locator.FindClosestPoint(p1)
    id2 = locator.FindClosestPoint(p2)
    if id1 < id2:
        return 1, 2
    else:
        return 2, 1

#!/usr/bin/env python

from common import *
import vtk
from scipy.interpolate import splrep, splev
import sys
import math
from vmtk import vtkvmtk
import numpy as np

def CreateParentArteryPatches(parentCenterlines, clipPoints, clipped=False):
    numberOfDaughterPatches = parentCenterlines.GetNumberOfCells()
    clipIds, numberOfPatchedCenterlinesPoints = ExtractPatchesIds(parentCenterlines, clipPoints, clipped)
    pnt = []

    patchedCenterlines = vtk.vtkPolyData()
    patchedCenterlinesPoints = vtk.vtkPoints()
    patchedCenterlinesCellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, numberOfPatchedCenterlinesPoints)

    numberOfCommonPatch = clipIds[0]+1
    patchedCenterlinesCellArray.InsertNextCell(numberOfCommonPatch)

    count = 0
    line = ExtractSingleLine(parentCenterlines, 0)
    getData = line.GetPointData().GetArray(radiusArrayName).GetTuple1
    for i in range(0, numberOfCommonPatch):
        patchedCenterlinesPoints.InsertNextPoint(line.GetPoint(i))
        patchedCenterlinesCellArray.InsertCellPoint(i)
        radiusArray.SetTuple1(i, getData(i))
        count += 1

    for j in range(numberOfDaughterPatches):
        cell = ExtractSingleLine(parentCenterlines, j)

        getData = cell.GetPointData().GetArray(radiusArrayName).GetTuple1
        numberOfCellPoints = cell.GetNumberOfPoints()
        #startId = clipIds[j]
        startId = clipIds[j+1]


        patchNumberOfPoints = numberOfCellPoints-startId
        patchedCenterlinesCellArray.InsertNextCell(patchNumberOfPoints)

        for i in range(startId, cell.GetNumberOfPoints()):
            point = cell.GetPoint(i)
            patchedCenterlinesPoints.InsertNextPoint(point)
            patchedCenterlinesCellArray.InsertCellPoint(count)
            radiusArray.SetTuple1(count, getData(i))
            count+=1

    patchedCenterlines.SetPoints(patchedCenterlinesPoints)
    patchedCenterlines.SetLines(patchedCenterlinesCellArray)
    patchedCenterlines.GetPointData().AddArray(radiusArray)

    return patchedCenterlines


def ExtractPatchesIds(parentCl, clipPts, clipped=False):
    clipIds = []
    numberOfPoints = 0
    N = clipPts.GetNumberOfPoints()

    upstreamPoint =  clipPts.GetPoint(0)
    downstreamPoint = clipPts.GetPoint(1)

    for j in range(parentCl.GetNumberOfCells()):
        cellLine = ExtractSingleLine(parentCl, j)
        locator = get_locator(cellLine)

        upId = locator.FindClosestPoint(upstreamPoint)
        downId = locator.FindClosestPoint(downstreamPoint)
        
        if j == 0:
            if clipped:
                clipIds.append(upId-1)
                clipIds.append(downId+1)
            else:
                clipIds.append(upId)
                clipIds.append(downId)
            numberOfPoints += upId + 1
            numberOfPoints += cellLine.GetNumberOfPoints() - downId 
        else:
            if clipped:
                clipIds.append(downId+1)
            else:
                clipIds.append(downId)
            numberOfPoints += cellLine.GetNumberOfPoints() - downId  
 
    return clipIds, numberOfPoints


def InterpolatePatchCenterlines(patchCenterlines, parentCenterlines,
                                additionalPoints, additionalPointIds, version_):
#def InterpolatePatchCenterlines(patchCenterlines, parentCenterlines,
#                                additionalPoint, version_):
    
    if additionalPoints is not None:
        additionalPointIds = []
        for i in range(parentCenterlines.GetNumberOfCells()):
            additionalPointIds.append([])
            for j in range(len(additionalPoints)):
                line = ExtractSingleLine(parentCenterlines, i)
                additionalPointIds[i].append(line.FindPoint(additionalPoints[j]))
    else:
        additionalPointIds = ["" for i in range(parentCenterlines.GetNumberOfCells())]
    
    interpolatedLines = vtk.vtkPolyData()
    interpolatedPoints = vtk.vtkPoints()
    interpolatedCellArray = vtk.vtkCellArray()
    radiusArray = get_vtk_array(radiusArrayName, 1, parentCenterlines.GetNumberOfCells())

    pointsInserted = 0
    interpolatedCellArray.Initialize()

    for i in range(parentCenterlines.GetNumberOfCells()):
        startingCell = vtk.vtkGenericCell()
        endingCell = vtk.vtkGenericCell()

        numberOfInterpolationPoints = parentCenterlines.GetCell(i).GetNumberOfPoints()

        patchCenterlines.GetCell(0, startingCell)
        patchCenterlines.GetCell(i+1, endingCell)

        end_points = []
        start_points = []

        for j in range(startingCell.GetNumberOfPoints()):
            start_points.append(startingCell.GetPoints().GetPoint(j))

        for k in range(endingCell.GetNumberOfPoints()):
            end_points.append(endingCell.GetPoints().GetPoint(k))

        #import matplotlib.pyplot as plt
        #from mpl_toolkits.mplot3d import Axes3D
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')

        end_p = np.asarray(end_points)
        start_p = np.asarray(start_points)

        #leg = str(version_) if lower is None else str(version_) + " " + lower
        #ax.plot(end_p[:,0], end_p[:,1], end_p[:, 2], label="end" + leg)
        #ax.plot(start_p[:,0], start_p[:,1], start_p[:,2], label="start" + leg)
        #if additionalPoint is not None:
        #    ax.plot([additionalPoint[0]], [additionalPoint[1]], [additionalPoint[2]], "o")
        #plt.legend()
        #plt.show()

        if version_:
            splinePoints = InterpolateSpline(startingCell, endingCell, additionalPoints)
        else:
            splinePoints = InterpolateTwoCells(startingCell, endingCell, \
                                               numberOfInterpolationPoints, \
                                               additionalPointIds,
                                               additionalPoints)
                                               #additionalPointIds[i],
                                               #additionalPoints)

        interpolatedCellArray.InsertNextCell(splinePoints.GetNumberOfPoints())
        for j in range(splinePoints.GetNumberOfPoints()):
            interpolatedPoints.InsertNextPoint(splinePoints.GetPoint(j))
            interpolatedCellArray.InsertCellPoint(pointsInserted + j)
        pointsInserted += splinePoints.GetNumberOfPoints()

    interpolatedLines.SetPoints(interpolatedPoints)
    interpolatedLines.SetLines(interpolatedCellArray)

    attributeFilter = vtkvmtk.vtkvmtkCenterlineAttributesFilter()
    if version < 6:
        attributeFilter.SetInput(interpolatedLines)
    else:
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
    N = 100#number of knots
    num_centerline_points = 3

    for i in range(num_centerline_points -1, -1, -1):
        points.append(get_startCell.GetPoint(num_start - n*i - 1))

    if additionalPoint is not None:
        for i in range(len(additionalPoint)):
            points.append(additionalPoint[i])

    for i in range(num_centerline_points):
        points.append(get_endCell.GetPoint(i*n))

    #from IPython import embed; embed()
    curv_coor = np.zeros(len(points))
    for i in range(len(points)-1):
        curv_coor[i+1] = curv_coor[i] + math.sqrt(distance(points[i], points[i+1]))

    # FIXME: Test the effect of adding an extra distance
    # Increse the longest distances
    #base_line = curv_coor[0]
    #for i, c in enumerate(curv_coor):
    #    if c > 5*base_line:
    #        print("Making distance longer, prev:", c, "  new:", 1.2*c)
    #        curv_coor[i] = c*1.2

    points = np.asarray(points)

    fx = splrep(curv_coor, points[:,0], k=3)
    fy = splrep(curv_coor, points[:,1], k=3)
    fz = splrep(curv_coor, points[:,2], k=3)

    curv_coor = np.linspace(curv_coor[0], curv_coor[-1], N)
    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    tmp = []
    for i in range(num_start - n*num_centerline_points):
        tmp.append(get_startCell.GetPoint(i))

    for j in range(N):
        tmp.append([fx_[j], fy_[j], fz_[j]])

    for k in range(n*num_centerline_points, num_end):
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
        for i in range(len(additionalPoint)):
            xspline.AddPoint(float(additionalPointId[i]), additionalPoint[i][0])
            yspline.AddPoint(float(additionalPointId[i]), additionalPoint[i][1])
            zspline.AddPoint(float(additionalPointId[i]), additionalPoint[i][2])

    for i in range(numberOfEndCellPoints):
        point = endCell.GetPoints().GetPoint(i)
        index = float(endCellFirstId + i)
        xspline.AddPoint(index, point[0])
        yspline.AddPoint(index, point[1])
        zspline.AddPoint(index, point[2])

    xspline.Compute()
    yspline.Compute()
    zspline.Compute()

    # Penelty is to avoid an adbrupt transition between interpolated and parent
    # centerline
    points.SetNumberOfPoints(numberOfSplinePoints - 1)
    penelty = 0
    for i in range(numberOfSplinePoints):
        if i == numberOfStartCellPoints:
            penelty = 1
            continue
        points.SetPoint(i - penelty, xspline.Evaluate(float(i)), yspline.Evaluate(float(i)), zspline.Evaluate(float(i)))

    return points

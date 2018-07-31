#!/usr/bin/env python

## Program:   AneuTool
## Module:    ToolRepairSTL.py
## Language:  Python
## Date:      $Date: 2016/17/04 00:00:00 $
## Version:   $Revision: 0.0.1 $
## Author:    Christophe Chnafa

##   Copyright (c) Christophe Chnafa. All rights reserved.

import argparse
import vtk
import math

from ImportData import loadFile


def isNaN(num):
    return num != num


def computeQualityMesh(mesh, measure):
    quality = vtk.vtkMeshQuality()
    quality.SetInput(mesh)
    #quality.SetTetQualityMeasure(measure) # TODO, test if 3d 
    quality.SetTriangleQualityMeasureToAspectRatio()
    quality.SetTriangleQualityMeasureToArea()
    quality.Update()
    return quality


def DumpQualityStats(iq, arrayname):
    an = iq.GetOutput().GetFieldData().GetArray(arrayname)
    cardinality = an.GetComponent(0, 4)
    range = list()
    range.append(an.GetComponent(0, 0))
    range.append(an.GetComponent(0, 2))
    average = an.GetComponent(0, 1)
    stdDev = math.sqrt(math.fabs(an.GetComponent(0, 3)))
    outStr = '%s%g%s%g%s%g\n%s%g%s%g' % (
            '  cardinality: ', cardinality,
            '  , range: ', range[0], '  -  ', range[1],
            '  average: ', average, '  , standard deviation: ', stdDev)
    print(outStr)
    return outStr


def surfaceOverview(mesh):
    nTriangles = mesh.GetNumberOfCells()
    nPoints = mesh.GetNumberOfPoints()
    print("> --- Surface overview:")
    print(("> Total number of triangles: %s." % nTriangles))
    print(("> Total number of points: %s." % nPoints))
    print(">")


def foundAndDeleteNaNTriangles(mesh):
    ctrNaN = 0
    foundNaN = False
    nTriangles = mesh.GetNumberOfCells()
    print("> --- Check the surface.")
    # The links from points to cells need to be build.
    mesh.BuildLinks()
    for i in range(0, nTriangles):
        killThisTriangle = False
        nPointsForThisCell = mesh.GetCell(i).GetPoints().GetNumberOfPoints()
        if  nPointsForThisCell> 3:
            print("> WARNING: found Cell with more than 3 points: there is more than triangles.")
        for j in range(0, nPointsForThisCell):
            x = [0.0, 0.0, 0.0]
            mesh.GetCell(i).GetPoints().GetPoint(j, x)
            if isNaN(x[0]) | isNaN(x[1]) | isNaN(x[2]):
                ctrNaN += 1
                killThisTriangle = True
        if killThisTriangle == True:
            mesh.DeleteCell(i)
    mesh.RemoveDeletedCells()
    print(("> Found %s NaN cells." % ctrNaN))
    print(">")
    if ctrNaN > 0:
        foundNaN = True

    return(foundNaN)


def cleanTheSurface(mesh):
    print("> --- Cleaning the surface.")
    cleanPolyData = vtk.vtkCleanPolyData()
    if vtk.VTK_MAJOR_VERSION <= 5:
        cleanPolyData.SetInput(mesh)
    else:
        cleanPolyData.SetInputData(mesh)
    #cleanPolyData.PointMergingOn()
    cleanPolyData.PointMergingOff() #point locator will not be used, and points 
                                    #that are not used by any cells will be eliminated, 
                                    #but never merged.
    #cleanPolyData.PieceInvariantOn()
    #cleanPolyData.ConvertStripsToPolysOff()
    # In VTK the tolerance is defined as a fraction 
    # of the bounding box length.
    tol = 0.0 #0.0005
    cleanPolyData.SetTolerance(tol)
    cleanPolyData.Update()

    cleanPolyData.Update()
    outputPolyData = cleanPolyData.GetOutput()
    
    print("> Done.")
    print("> ")

    return(outputPolyData)


def closeAllTheHolesOnTheSurface(mesh):
    print("> --- Closing the surface.")
    fillHolesFilter = vtk.vtkFillHolesFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
        fillHolesFilter.SetInput(mesh)
    else:
        fillHolesFilter.SetInputConnection(mesh.GetOutputPort())
    fillHolesFilter.SetHoleSize(100000000000)
    fillHolesFilter.Update()

    outputPolyData = fillHolesFilter.GetOutput()
    #outputPolyData.Update()
    print("> Done.")
    print("> ")

    return(outputPolyData)


def cleanTheSurfaceUnstructGrid(mesh):
    print("> --- Cleaning the surface.")
    cleanPolyData = vtk.vtkCleanUnstructuredGrid()
    if vtk.VTK_MAJOR_VERSION <= 5:
        cleanPolyData.SetInput(mesh)
    else:
        cleanPolyData.SetInputConnection(mesh.GetOutputPort())
    cleanPolyData.Update()
    outputPolyData = cleanPolyData.GetOutput()
    outputPolyData.Update()
    print("> Done.")
    print("> ")

    return(outputPolyData)


def writeSTLfile(outputPolyData, outputFileName, ascii):
    print("> --- Write the clean STL file as:", outputFileName)
    print(">")
    stlWriter = vtk.vtkSTLWriter()
    stlWriter.SetFileName(outputFileName)
    if not(ascii):
        stlWriter.SetFileTypeToBinary()
    if vtk.VTK_MAJOR_VERSION <= 5:
        stlWriter.SetInput(outputPolyData)
    else:
        stlWriter.SetInputConnection(outputPolyData.GetOutputPort())
    stlWriter.Write()


def checkIfThereIsHoles(outputBisPolyData):
    featureEdges = vtk.vtkFeatureEdges()
    featureEdges.FeatureEdgesOff()
    featureEdges.BoundaryEdgesOn()
    featureEdges.NonManifoldEdgesOn()
    featureEdges.SetInput(outputBisPolyData)
    featureEdges.Update()
    numberOfOpenEdges = featureEdges.GetOutput().GetNumberOfCells()
    if numberOfOpenEdges > 0:
        print("> The surface is open. It might be only the outlets if you did not try to close the surface.")
    else:
        print("> The surface is closed, there is no holes.")
    print(">")

def checkIfThereIsNonTriangleCells(outputBisPolyData):
    ctr = 0
    for i in range(0, outputBisPolyData.GetNumberOfCells()):
        if outputBisPolyData.GetCellType(i) != 5:
            #print outputBisPolyData.GetCellType(i)
            ctr += 1
    if ctr > 0:
        print("> ERROR: There are elements that are not triangles")


def repairSTL(inputFileName, outputFileName, closeHoles, checkHoles, checkQuality, ascii):
    # Open file.
    mesh = loadFile(inputFileName)

    if checkQuality == True:
        quality = computeQualityMesh(mesh, vtk.VTK_QUALITY_ASPECT_RATIO)
        #q = qualityFilter.GetOutput().GetCellData().GetArray("Quality")

        print(DumpQualityStats(quality, 'Quality'))
        # quality_min = quality.GetOutput().GetFieldData()\
        # .GetArray('Mesh Tetrahedron Quality').GetComponent(0, 0)
        # quality_max = quality.GetOutput().GetFieldData()\
        # .GetArray('Mesh Tetrahedron Quality').GetComponent(0, 2)
        # print quality_min
        # print quality_max

    # Check the number of pts and cells.
    surfaceOverview(mesh)

    # Check if there is NaN coordinates and clean the surface.
    foundAndDeleteNaNTriangles(mesh)

    # Check the number of pts and cells.
    surfaceOverview(mesh)

    # Clean the surface.
    outputPolyData = cleanTheSurface(mesh)

    # Check the number of pts and cells.
    surfaceOverview(outputPolyData)

    # Re check the surface.
    foundNaN = foundAndDeleteNaNTriangles(outputPolyData)


    if foundNaN == True:
        print("> ERROR: There is still an issue that cannot be fixed with this tool.")
    else:
        if closeHoles == True:
            outputBisPolyData = closeAllTheHolesOnTheSurface(outputPolyData)
            surfaceOverview(outputBisPolyData)
        else:
            outputBisPolyData = outputPolyData            

        if checkHoles == True:
            checkIfThereIsHoles(outputBisPolyData)

        checkIfThereIsNonTriangleCells(outputBisPolyData)

    if not(outputFileName == ''):
        writeSTLfile(outputBisPolyData, outputFileName, ascii)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = "ToolRepairSTL: correct an STL surface with NaN points and/or redundant nodes.")
    parser.add_argument('-i', '--input', type = str, required = True, 
        dest = 'inputFileName',
        help = "Input file containing the surface data in a STL format.")
    parser.add_argument('-o', '--output', type = str, required = False, 
        default = '', dest = 'outputFileName',
        help = "Output file containing the cleaned surface data in a STL format.")
    parser.add_argument('-p', '--patchesholes', required = False, 
        default = False, 
        dest='closeHoles', action = "store_true", 
        help = "Set this argument to true if you want to patch up the potential holes of the surface.")
    parser.add_argument('-c', '--checkholes', required = False, 
        default = False, 
        dest='checkholes', action = "store_true", 
        help = "Set this argument to true if you want to check up the potential holes of the surface.")
    parser.add_argument('-a', '--ascii', required = False, 
        default = False, 
        dest='ascii', action = "store_true", 
        help = "Set this argument to true if you want the output STL wrote in ASCII.")
    parser.add_argument('-q', '--checkquality', required = False, 
        default = False, 
        dest='checkQuality', action = "store_true", 
        help = "Set this argument to true if you want to check up the potential holes of the surface.")
    

    args = parser.parse_args()
    repairSTL(args.inputFileName, args.outputFileName, args.closeHoles, args.checkholes, args.checkQuality, args.ascii)

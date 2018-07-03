#!/usr/binY/env python

## Program:   VMTK
## Module:    $RCSfile: vmtkpointselector.py,v $
## Language:  Python
## Date:      $Date: 2006/07/17 09:52:56 $
## Version:   $Revision: 1.20 $

##   Copyright (c) Luca Antiga, David Steinman. All rights reserved.
##   See LICENCE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.


import vtk
import sys
from os import path

from vmtk import vtkvmtk
from vmtk import vmtkrenderer
from vmtk import pypes

vmtkpointselector = 'vmtkPointselector'


class vmtkSeedSelector():

    def __init__(self):
        self._Surface = None
        self._SeedIds = None
        self._SourceSeedIds = vtk.vtkIdList()
        self._TargetSeedIds = vtk.vtkIdList()
        self.PrintError = None
        self.PrintLog = None
        self.InputText = None
        self.OutputText = None
        self.InputInfo = None

    def SetSurface(self,surface):
        self._Surface = surface

    def GetSurface(self):
        return self._Surface

    def GetSourceSeedIds(self):
        return self._SourceSeedIds

    def GetTargetSeedIds(self):
        return self._TargetSeedIds

    def Execute(self):
        pass


class vmtkPickPointSeedSelector(vmtkSeedSelector):

    def __init__(self):
        vmtkSeedSelector.__init__(self)
        self.PickedSeedIds = vtk.vtkIdList()
        self.PickedSeeds = vtk.vtkPolyData()
        self.vmtkRenderer = None
        self.OwnRenderer = 0
        self.Script = None

    def UndoCallback(self, obj):
        self.InitializeSeeds()
        self.PickedSeeds.Modified()
        self.vmtkRenderer.RenderWindow.Render()

    def PickCallback(self, obj):
        picker = vtk.vtkCellPicker()
        picker.SetTolerance(1E-4 * self._Surface.GetLength())
        eventPosition = self.vmtkRenderer.RenderWindowInteractor.GetEventPosition()
        result = picker.Pick(float(eventPosition[0]),float(eventPosition[1]),0.0,self.vmtkRenderer.Renderer)
        if result == 0:
            return
        pickPosition = picker.GetPickPosition()
        pickedCellPointIds = self._Surface.GetCell(picker.GetCellId()).GetPointIds()
        minDistance = 1E10
        pickedSeedId = -1
        for i in range(pickedCellPointIds.GetNumberOfIds()):
            distance = vtk.vtkMath.Distance2BetweenPoints(pickPosition,self._Surface.GetPoint(pickedCellPointIds.GetId(i)))
            if distance < minDistance:
                minDistance = distance
                pickedSeedId = pickedCellPointIds.GetId(i)
        if pickedSeedId == -1:
            pickedSeedId = pickedCellPointIds.GetId(0)
        self.PickedSeedIds.InsertNextId(pickedSeedId)
        point = self._Surface.GetPoint(pickedSeedId)
        self.PickedSeeds.GetPoints().InsertNextPoint(point)
        self.PickedSeeds.Modified()
        self.vmtkRenderer.RenderWindow.Render()

    def InitializeSeeds(self):
        self.PickedSeedIds.Initialize()
        self.PickedSeeds.Initialize()
        seedPoints = vtk.vtkPoints()
        self.PickedSeeds.SetPoints(seedPoints)

    def Execute(self):

        if (self._Surface == None):
            self.PrintError('vmtkPickPointSeedSelector Error: Surface not set.')
            return

        self._SourceSeedIds.Initialize()
        self._TargetSeedIds.Initialize()

        if not self.vmtkRenderer:
            self.vmtkRenderer = vmtkrenderer.vmtkRenderer()
            self.vmtkRenderer.Initialize()
            self.OwnRenderer = 1

        self.vmtkRenderer.RegisterScript(self.Script) 

        glyphs = vtk.vtkGlyph3D()
        glyphSource = vtk.vtkSphereSource()
        glyphs.SetInput(self.PickedSeeds)
        glyphs.SetSource(glyphSource.GetOutput())
        glyphs.SetScaleModeToDataScalingOff()
        glyphs.SetScaleFactor(self._Surface.GetLength()*0.01)
        glyphMapper = vtk.vtkPolyDataMapper()
        glyphMapper.SetInput(glyphs.GetOutput())
        self.SeedActor = vtk.vtkActor()
        self.SeedActor.SetMapper(glyphMapper)
        self.SeedActor.GetProperty().SetColor(1.0,0.0,0.0)
        self.SeedActor.PickableOff()
        self.vmtkRenderer.Renderer.AddActor(self.SeedActor)

        ##self.vmtkRenderer.RenderWindowInteractor.AddObserver("KeyPressEvent", self.KeyPressed)
        self.vmtkRenderer.AddKeyBinding('u','Undo.',self.UndoCallback)
        self.vmtkRenderer.AddKeyBinding('space','Add points.',self.PickCallback)
        surfaceMapper = vtk.vtkPolyDataMapper()
        surfaceMapper.SetInput(self._Surface)
        surfaceMapper.ScalarVisibilityOff()
        surfaceActor = vtk.vtkActor()
        surfaceActor.SetMapper(surfaceMapper)
        surfaceActor.GetProperty().SetOpacity(1.0)

        self.vmtkRenderer.Renderer.AddActor(surfaceActor)

        self.InputInfo('Please position the mouse and press space to add source points, \'u\' to undo\n')

        any = 0
        while any == 0:
            self.InitializeSeeds()
            self.vmtkRenderer.Render()
            any = self.PickedSeedIds.GetNumberOfIds()
        self._SourceSeedIds.DeepCopy(self.PickedSeedIds)

        self.InputInfo('Please position the mouse and press space to add target points, \'u\' to undo\n')

        any = 0
        while any == 0:
            self.InitializeSeeds()
            self.vmtkRenderer.Render()
            any = self.PickedSeedIds.GetNumberOfIds()
        self._TargetSeedIds.DeepCopy(self.PickedSeedIds)

        if self.OwnRenderer:
            self.vmtkRenderer.Deallocate()


class vmtkNonManifoldSurfaceChecker(object):

    def __init__(self):

        self.Surface = 0
        
        self.NumberOfNonManifoldEdges = 0
        self.Report = 0
        self.NonManifoldEdgePointIds = vtk.vtkIdList()

        self.PrintError = None

    def Execute(self):

        if (self.Surface == 0):
            self.PrintError('NonManifoldSurfaceChecker error: Surface not set')
            return

        self.NonManifoldEdgesFound = 0
        self.Report = ''
        self.NonManifoldEdgePointIds.Initialize()
        
        neighborhoods = vtkvmtk.vtkvmtkNeighborhoods()
        neighborhoods.SetNeighborhoodTypeToPolyDataManifoldNeighborhood()
        neighborhoods.SetDataSet(self.Surface)
        neighborhoods.Build()

        neighborCellIds = vtk.vtkIdList()
        cellPointIds = vtk.vtkIdList()

        self.Surface.BuildCells()
        self.Surface.BuildLinks(0)
        self.Surface.Update()

        numberOfNonManifoldEdges = 0

        for i in range(neighborhoods.GetNumberOfNeighborhoods()):

            neighborhood = neighborhoods.GetNeighborhood(i)
            
            for j in range(neighborhood.GetNumberOfPoints()):
                
                neighborId = neighborhood.GetPointId(j)
                
                if (i<neighborId):
                    
                    neighborCellIds.Initialize()
                    self.Surface.GetCellEdgeNeighbors(-1,i,neighborId,neighborCellIds)
                    
                    if (neighborCellIds.GetNumberOfIds()>2):

                        numberOfNonManifoldEdges = numberOfNonManifoldEdges + 1
                        
                        self.Report = self.Report +  "Non-manifold edge found" + str(i) + ' ' + str(neighborId) + '.\n'

                        self.NonManifoldEdgePointIds.InsertNextId(i)
                        self.NonManifoldEdgePointIds.InsertNextId(neighborId)


class vmtkPointselector(pypes.pypeScript):

    def __init__(self):

        pypes.pypeScript.__init__(self)
        
        self.Surface = None
        self.Output = None
        self.Centerlines = None
        self.SeedSelector = None
        self.SeedSelectorName = 'pickpoint'
        self.FlipNormals = 0
        self.CapDisplacement = 0.0
        self.RadiusArrayName = 'MaximumInscribedSphereRadius'
        self.CostFunction = '1/R'
        self.AppendEndPoints = 0
        self.CheckNonManifold = 0
        
        self.Resampling = 0
        self.ResamplingStepLength = 1.0
        self.SimplifyVoronoi = 0

        self.EikonalSolutionArrayName = 'EikonalSolution'
        self.EdgeArrayName = 'EdgeArray'
        self.EdgePCoordArrayName = 'EdgePCoordArray'
        self.CostFunctionArrayName = 'CostFunctionArray'

        self.UseTetGen = 0
        self.TetGenDetectInter = 1

        self.DelaunayTessellation = None
        self.VoronoiDiagram = None
        self.PoleIds = None

        self.DelaunayTolerance = 0.001

        self.SourceIds = []
        self.TargetIds = []
        self.SourcePoints = []
        self.TargetPoints = []

        self.vmtkRenderer = None
        self.OwnRenderer = 0

        self.SetScriptName('vmtkpointselector')
        self.SetScriptDoc('Choose points on a geometry that will be stored for later use')
        self.SetInputMembers([['Surface', 'i', 'vtkPolyData', 1, '',
                               'the input surface', 'vmtksurfacereader'],
                               ["Output", 'opath', 'str', 1, "", "Filepath to store points"]])
        self.SetOutputMembers([['Centerlines','o','vtkPolyData',1,'','the output centerlines',
                                'vmtksurfacewriter']])
    
    def PrintProgress(self,obj,event):
        self.OutputProgress(obj.GetProgress(),10)

    def Execute(self):
        self.SeedSelectorName = "pickpoint"

        if self.Surface == None:
            self.PrintError('Error: No input surface.')
        
        if self.CheckNonManifold:
            self.PrintLog('NonManifold check.')
            nonManifoldChecker = vmtkNonManifoldSurfaceChecker()
            nonManifoldChecker.Surface = self.Surface
            nonManifoldChecker.PrintError = self.PrintError
            nonManifoldChecker.Execute()

            if (nonManifoldChecker.NumberOfNonManifoldEdges > 0):
                self.PrintLog(nonManifoldChecker.Report)
                return

        self.vmtkRenderer = vmtkrenderer.vmtkRenderer()
        self.vmtkRenderer.Initialize()
        self.OwnRenderer = 1

        self.PrintLog('Cleaning surface.')
        surfaceCleaner = vtk.vtkCleanPolyData()
        surfaceCleaner.SetInput(self.Surface)
        surfaceCleaner.Update()

        self.PrintLog('Triangulating surface.')
        surfaceTriangulator = vtk.vtkTriangleFilter()
        surfaceTriangulator.SetInput(surfaceCleaner.GetOutput())
        surfaceTriangulator.PassLinesOff()
        surfaceTriangulator.PassVertsOff()
        surfaceTriangulator.Update()

        self.PrintLog('Capping surface.')
        surfaceCapper = vtkvmtk.vtkvmtkCapPolyData()
        surfaceCapper.SetInput(surfaceTriangulator.GetOutput())
        surfaceCapper.SetDisplacement(self.CapDisplacement)
        surfaceCapper.SetInPlaneDisplacement(self.CapDisplacement)
        surfaceCapper.Update()
        centerlineInputSurface = surfaceCapper.GetOutput()
        capCenterIds = surfaceCapper.GetCapCenterIds()

        self.SeedSelector = vmtkPickPointSeedSelector()
        self.SeedSelector.vmtkRenderer = self.vmtkRenderer
        self.SeedSelector.Script = self

        self.SeedSelector.SetSurface(centerlineInputSurface)
        self.SeedSelector.InputInfo = self.InputInfo
        self.SeedSelector.InputText = self.InputText
        self.SeedSelector.OutputText = self.OutputText
        self.SeedSelector.PrintError = self.PrintError
        self.SeedSelector.PrintLog = self.PrintLog
        self.SeedSelector.Execute()

        inletSeedIds = self.SeedSelector.GetSourceSeedIds()
        outletSeedIds = self.SeedSelector.GetTargetSeedIds()

        point_source = []
        points_target = []
        s = ""
        text_target = "target%s: %s\n"
        text_source = "source: %s\n"
        
        for i in range(inletSeedIds.GetNumberOfIds()):
            point_source.append(centerlineInputSurface.GetPoint(inletSeedIds.GetId(i)))
            s += text_source % (str(point_source[-1]))

        for i in range(outletSeedIds.GetNumberOfIds()):
            points_target.append(centerlineInputSurface.GetPoint(outletSeedIds.GetId(i)))
            s += text_target % (i, points_target[-1])

        f = open(path.join(self.Output, "points"), "w")
        f.write(s)
        f.close()


if __name__=='__main__':
    main = pypes.pypeMain()
    main.Arguments = sys.argv
    main.Execute()

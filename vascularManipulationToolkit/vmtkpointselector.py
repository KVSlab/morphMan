#!/usr/bin/env python

import vtk
import sys
from os import path
from vmtk.vmtkscripts import vmtkRenderer
from vmtk import vmtkrenderer

version = vtk.vtkVersion().GetVTKMajorVersion()

class VtkText:
    def __init__(self, guiText=""):
        self.text = vtk.vtkTextActor()
        self.text.SetInput(guiText)
        textProperties = self.text.GetTextProperty()
        #textProperties.SetFontFamilyToArial()
        textProperties.SetFontSize(15)
        textProperties.SetColor(1, 1, 1)
        self.text.SetDisplayPosition(10, 2)


class vmtkSeedSelector():
    def __init__(self):
        self._Surface = None
        self._SeedIds = None
        self._TargetSeedIds = vtk.vtkIdList()

    def SetSurface(self, surface):
        self._Surface = surface

    def GetSurface(self):
        return self._Surface

    def GetTargetSeedIds(self):
        return self._TargetSeedIds

    def Execute(self, mode=None):
        pass


class vmtkPickPointSeedSelector(vmtkSeedSelector):
    def __init__(self):
        vmtkSeedSelector.__init__(self)
        self.PickedSeedIds = vtk.vtkIdList()
        self.PickedSeeds = vtk.vtkPolyData()
        self.vmtkRenderer = None
        self.OwnRenderer = 0
        self.Script = None
        self.Mode = None

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

        self._TargetSeedIds.Initialize()

        if not self.vmtkRenderer:
            self.vmtkRenderer = vmtkrenderer.vmtkRenderer()
            self.vmtkRenderer.Initialize()
            self.OwnRenderer = 1

        #self.vmtkRenderer.RegisterScript(self.Script) 

        glyphs = vtk.vtkGlyph3D()
        glyphSource = vtk.vtkSphereSource()
        if version < 6:
            glyphs.SetInput(self.PickedSeeds)
            glyphs.SetSource(glyphSource.GetOutput())
        else:
            glyphs.SetInputData(self.PickedSeeds)
            #glyphs.SetSourceData(glyphSource.GetOutput())
            glyphs.SetSourceConnection(glyphSource.GetOutputPort())

        glyphs.SetScaleModeToDataScalingOff()
        glyphs.SetScaleFactor(self._Surface.GetLength()*0.01)
        glyphMapper = vtk.vtkPolyDataMapper()
        if version < 6:
            glyphMapper.SetInput(glyphs.GetOutput())
        else:
            # NOTE: This fixed the red spheres not showing up
            glyphMapper.SetInputData(glyphs.GetOutput())
        glyphMapper.SetInputConnection(glyphs.GetOutputPort())
        self.SeedActor = vtk.vtkActor()
        self.SeedActor.SetMapper(glyphMapper)
        self.SeedActor.GetProperty().SetColor(1.0, 0.0, 0.0)
        self.SeedActor.PickableOff()
        self.vmtkRenderer.Renderer.AddActor(self.SeedActor)

        self.vmtkRenderer.AddKeyBinding('u', 'Undo.', self.UndoCallback)
        self.vmtkRenderer.AddKeyBinding('space', 'Add points.', self.PickCallback)
        surfaceMapper = vtk.vtkPolyDataMapper()
        if version < 6:
            surfaceMapper.SetInput(self._Surface)
        else:
            surfaceMapper.SetInputData(self._Surface)
        surfaceMapper.ScalarVisibilityOff()
        surfaceActor = vtk.vtkActor()
        surfaceActor.SetMapper(surfaceMapper)
        surfaceActor.GetProperty().SetOpacity(1.0)

        self.vmtkRenderer.Renderer.AddActor(surfaceActor)
        if self.Mode == "stenosis":
            text = 'Place one point where stenosis is centered OR place two points to heal stenosis, \'u\' to undo\n'
        elif self.Mode == "outlet":
            text = 'Please select the two relevant outlets, \'u\' to undo\n'
        elif self.Mode == "aneurysm": 
            text = 'Please position the mouse and press space to add the top of the aneurysm, \'u\' to undo\n'
        else:
            text = "Place points, \'u\' to undo\n"

        guiText = VtkText(text)
        self.vmtkRenderer.Renderer.AddActor(guiText.text)

        any = 0
        while any == 0:
            self.InitializeSeeds()
            self.vmtkRenderer.Render()
            any = self.PickedSeedIds.GetNumberOfIds()
        self._TargetSeedIds.DeepCopy(self.PickedSeedIds)

        if self.OwnRenderer:
            self.vmtkRenderer.Deallocate()

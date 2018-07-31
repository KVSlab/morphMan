#!/usr/bin/env python

## Program:   AneuTools
## Module:    ImportData.py
## Language:  Python
## Date:      $Date: 2016/17/04 00:00:00 $
## Version:   $Revision: 0.0.1 $
## Author:    Christophe Chnafa

##   Copyright (c) Christophe Chnafa. All rights reserved.

import math
import numpy as np
import vtk

# Array names used by VMTK.
BLANKINGARRAYNAME = 'Blanking'
GROUPIDSARRAYNAME = 'GroupIds'
RADIUSARRAYNAME = 'MaximumInscribedSphereRadius'
SECTIONARRAYNAME = 'CenterlineSectionArea'


def loadFile(fileName):
    '''Load the given file, and return a vtkPolyData object for it. '''
    fileType = fileName[-3:]
    if fileType == '':
        raise RuntimeError('The file does not have an extension')
    if fileType == 'stl':
        reader = vtk.vtkSTLReader()
        reader.MergingOn()
    elif fileType == 'vtk':
        reader = vtk.vtkPolyDataReader()
    elif fileType == 'vtp':
        reader = vtk.vtkXMLPolyDataReader()
    elif fileType == 'vtu':
        reader = vtk.vtkUnstructuredGridReader()
    else:
        raise RuntimeError('Unknown file type %s' % fileType)
    reader.SetFileName(fileName)
    reader.Update()
    polyData = reader.GetOutput()

    return(polyData)

def GetMaxGroupId(centerline):
    maxGroupId = 0
    groupIdsArray = centerline.GetCellData().GetArray(GROUPIDSARRAYNAME)
    for i in range(0, centerline.GetNumberOfCells()):
        groupId = groupIdsArray.GetValue(i)
        if groupId > maxGroupId:
            maxGroupId = groupId

    return(maxGroupId)

def GetBlankedGroupsIdList(centerline):
    blankedGroupdsArray = []
    idx = 0
    groupIdsArray = centerline.GetCellData().GetArray(GROUPIDSARRAYNAME)
    blankingArray = centerline.GetCellData().GetArray(BLANKINGARRAYNAME)
    for i in range(0, centerline.GetNumberOfCells()):
        groupId = groupIdsArray.GetValue(i)
        blanking = blankingArray.GetValue(i)
        if blanking == 1:
            blankedGroupdsArray.append((idx, groupId)) 
        idx += 1

    return(blankedGroupdsArray)

def GetRedundantBlankedIdList(centerline, blankedGroupsIdList):
    ''' Get the redundant blanked segments.

    Nominaly, the x0 and x1 coordinates should be the same. Thus, the
    computed distance should be exactly 0.0.
    Sometimes branches are just different by a few and SHOULD be merged.
    Visual inspection of the centerline with Paraview when the Warning 
    message pop up is advised.
    Sometimes, the problem is more concerning. It is problematic that 
    VMTK split, sometimes, the branches with a segment with a blanked 
    segment, a lil' something and another blanked segment. In case of 
    wierd branch splittings, try to lower the tol. It is cheating since 
    branches that should not be merge will be merged anyway.

    '''
    currentX0 = [0.0, 0.0, 0.0]
    currentX1 = [0.0, 0.0, 0.0]
    otherX0 = [0.0, 0.0, 0.0]
    otherX1 = [0.0, 0.0, 0.0]
    minLength, maxLength = ComputeGeometricTolerance(centerline)
    tol = maxLength*0.1 #10**(-5) 
    redundantBranchesIndex = []
    for currentBranch in blankedGroupsIdList:
        currentBranchIndex = currentBranch[0]
        if currentBranchIndex in redundantBranchesIndex:
            continue
        centerline.GetCell(currentBranchIndex).GetPoints().GetPoint(0, currentX0)
        npts = centerline.GetCell(currentBranchIndex).GetPoints().GetNumberOfPoints()
        centerline.GetCell(currentBranchIndex).GetPoints().GetPoint(npts - 1, currentX1)
        for otherBranch in blankedGroupsIdList:
            otherBranchIndex = otherBranch[0]
            if otherBranchIndex == currentBranchIndex:
                continue
            if otherBranchIndex in redundantBranchesIndex:
                continue
            centerline.GetCell(otherBranchIndex).GetPoints().GetPoint(0, otherX0)
            otherNpts = centerline.GetCell(otherBranchIndex).GetPoints().GetNumberOfPoints()
            centerline.GetCell(otherBranchIndex).GetPoints().GetPoint(otherNpts - 1, otherX1)
            if vtk.vtkMath.Distance2BetweenPoints(currentX0,otherX0) < tol and \
               vtk.vtkMath.Distance2BetweenPoints(currentX1,otherX1) < tol:
                redundantBranchesIndex.append(otherBranchIndex)
                if vtk.vtkMath.Distance2BetweenPoints(currentX0,otherX0) > minLength:
                    print('WARNING: POTENTIAL ISSUE DURING THE MERGING OF REDUNDANTS BLANKED SEGMENTS.')
                    print('         A distance between segments is suspicious.') 
                    print('         The blanked segments of VTK Cell Id ' + repr(currentBranchIndex))
                    print('         and, VTK Cell Id ' + repr(otherBranchIndex) + ' will be merged.') 
                    print('         The segments points should be identical.  ')
                    print('         Please check if this action was expected.')
                    print()

    return(redundantBranchesIndex)

def IsArrayDefined(centerline, arrayName):
    nPointDataArrays = centerline.GetPointData().GetNumberOfArrays()
    found = False
    for i in range(0, nPointDataArrays):
        if arrayName == centerline.GetPointData().GetArrayName(i):
            found = True
            break
    nCellDataArrays = centerline.GetCellData().GetNumberOfArrays()
    for i in range(0, nCellDataArrays):
        if arrayName == centerline.GetCellData().GetArrayName(i):
            found = True
            break

    return(found)

def ComputeGeometricTolerance(centerline):
    '''Return the min and max length for branches.

    This routine compute the delta x minimum and maximum in the network. 
    The minimum length is used for the connectivity. The maximum length
    is used for merging the blanked segments.

    '''
    groupIdsArray = centerline.GetCellData().GetArray(GROUPIDSARRAYNAME)
    numberOfCells = centerline.GetNumberOfCells()
    minLength = 10000.0
    maxLength = 0.0
    numberOfCells = centerline.GetNumberOfCells()
    for i in range(0, numberOfCells):
        groupId = groupIdsArray.GetValue(i)
        npts =  centerline.GetCell(i).GetPoints().GetNumberOfPoints()
        for k in range(0, npts - 1):
            point0 = [0.0, 0.0, 0.0]
            point1 = [0.0, 0.0, 0.0]
            centerline.GetCell(i).GetPoints().GetPoint(k, point0)
            centerline.GetCell(i).GetPoints().GetPoint(k + 1, point1)
            if math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1)) < minLength:
                if math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1)) == 0.0:
                    continue
                minLength = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))
            if math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))>maxLength:
                maxLength = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))

    return minLength, maxLength

def ComputeGroupLength(centerline, branchGroupId):
    '''Return the mean length for branches of branchGroupId.'''

    groupLength = 0.0
    lengthWeightSum = 0.0
    groupIdsArray = centerline.GetCellData().GetArray(GROUPIDSARRAYNAME)
    numberOfCells = centerline.GetNumberOfCells()
    for i in range(0, numberOfCells):
        groupId = groupIdsArray.GetValue(i)
        if groupId != branchGroupId:
            continue
        length = ComputeBranchLength(centerline, i)
        groupLength += length
        lengthWeightSum += 1.0
    groupLength /= lengthWeightSum

    return groupLength

def ComputeBranchLength(centerline, branchId):
    '''Return the length for the branch of index branchId.'''

    length = 0.0
    npts = centerline.GetCell(branchId).GetPoints().GetNumberOfPoints()
    for k in range(0, npts - 1):
        point0 = [0.0, 0.0, 0.0]
        point1 = [0.0, 0.0, 0.0]
        centerline.GetCell(branchId).GetPoints().GetPoint(k, point0)
        centerline.GetCell(branchId).GetPoints().GetPoint(k + 1, point1)
        length += math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))

    return length

def ComputeGroupRadius(centerline, branchGroupId):
    '''Return the mean radius of a group.

    The mean radius is computed using the hydraulic resistance for 
    the branches of branchGroupId. See Chnafa et al. International 
    Journal for Numerical Methods in Biomedical Engineering, 2016.

    '''
    groupRadius = 0.0
    radiusWeightSum = 0.0
    groupIdsArray = centerline.GetCellData().GetArray(GROUPIDSARRAYNAME)
    numberOfCells = centerline.GetNumberOfCells()
    for i in range(0, numberOfCells):
        if groupIdsArray.GetValue(i) != branchGroupId:
            continue
        groupRadius += ComputeBranchRadius(centerline, i)
        radiusWeightSum += 1.0
    groupRadius /= radiusWeightSum

    return groupRadius

def ComputeBranchRadius(centerline, branchId):
    '''Return the radius for a branch with index branchId. '''
    length = 0.0
    resistance = 0.0
    radiusArray = centerline.GetPointData().GetArray(RADIUSARRAYNAME)
    npts = centerline.GetCell(branchId).GetPoints().GetNumberOfPoints()
    for k in range(0, npts - 1):
        point0 = [0.0, 0.0, 0.0]
        point1 = [0.0, 0.0, 0.0]
        centerline.GetCell(branchId).GetPoints().GetPoint(k, point0)
        centerline.GetCell(branchId).GetPoints().GetPoint(k + 1, point1)
        dx = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))
        length += dx
        r = radiusArray.GetComponent(centerline.GetCell(branchId).GetPointId(k),0)
        resistance += dx / r**(4.0)
    branchRadius = (length / resistance)**(0.25)

    return branchRadius

def GetListsUniqueBlankedBranches(blankedGroupsIdList, redundantBlankedBranchesIdList):
    blankedGroupsIndex = []
    blankedUniqueBranchesIndex = []
    for i in range(0, len(blankedGroupsIdList)):
        if blankedGroupsIdList[i][0] in redundantBlankedBranchesIdList:
            continue
        blankedGroupsIndex.append(blankedGroupsIdList[i][1])
        blankedUniqueBranchesIndex.append(blankedGroupsIdList[i][0])

    return blankedGroupsIndex, blankedUniqueBranchesIndex

def SetNetworkStructure(centerline, network, verboseprint,
    isConnectivityNeeded = True, isRadiusInletNeeded = True):
    '''Fills a network structure with a vtkPolyData object.

    Each element has an unique index. The groups length and radius are
    averaged out. When the element is not not part of the bifurcation 
    reference system, it is possible that multiple end points exist. They are
    all saved. For the complet formalism of the bifurcation reference system 
    see Antiga, L., & Steinman, D. A. (2004). Robust and objective 
    decomposition and mapping of bifurcating vessels. Medical Imaging, 
    IEEE Transactions on, 23(6), 704-713.

    '''
    verboseprint("> Filling the network structure with the raw data.")
    if not(IsArrayDefined(centerline, GROUPIDSARRAYNAME)):
        from vmtk import vmtkscripts
        centerlines = vmtkscripts.vmtkBranchExtractor()
        centerlines.Centerlines = centerline
        centerlines.RadiusArrayName = RADIUSARRAYNAME
        centerlines.Execute()
        centerline = centerlines.Centerlines

    # Treat the splitted centerline.
    maxGroupId = GetMaxGroupId(centerline)
    blankedGroupsIdList = GetBlankedGroupsIdList(centerline)
    redundantBlankedBranchesIdList = GetRedundantBlankedIdList(centerline, 
        blankedGroupsIdList)
    blankedGroupsIndex, blankedUniqueBranchesIndex = \
        GetListsUniqueBlankedBranches(blankedGroupsIdList, 
        redundantBlankedBranchesIdList)
    indexUniqueBranches = 0
    for i in range(0, maxGroupId + 1):
        x0List = []
        x1List = []
        VtkCellIdList = []
        VtkGroupIdList = []
        if i in blankedGroupsIndex:
            for j in range(0,len(blankedGroupsIndex)):    
                if blankedGroupsIndex[j] == i:
                    el = Element(Id = indexUniqueBranches)
                    el.SetMeanRadius(ComputeBranchRadius(centerline,
                        blankedUniqueBranchesIndex[j]))
                    el.SetLength(ComputeBranchLength(centerline,
                        blankedUniqueBranchesIndex[j]))
                    el.SetBlanking(1)
                    npts = centerline.GetCell(blankedUniqueBranchesIndex[j]). \
                        GetPoints().GetNumberOfPoints()
                    x0 = [0.0, 0.0, 0.0]
                    x1 = [0.0, 0.0, 0.0]
                    centerline.GetCell(blankedUniqueBranchesIndex[j]). \
                        GetPoints().GetPoint(0, x0)
                    centerline.GetCell(blankedUniqueBranchesIndex[j]). \
                        GetPoints().GetPoint(npts - 1, x1)
                    x0List = []
                    x1List = []
                    VtkCellIdList =[]
                    VtkGroupIdList = []
                    x1List.append(x1)
                    x0List.append(x0)
                    el.SetInOutPointsCoordinates(x0List, x1List)
                    VtkCellIdList.append(blankedUniqueBranchesIndex[j])
                    VtkGroupIdList.append(i)
                    el.SetVtkGroupIdList(VtkGroupIdList)
                    el.SetVtkCellIdList(VtkCellIdList)
                    network.AddElement(el)
                    verboseprint("> Edge Id " + repr(indexUniqueBranches) 
                        + ", Length = " + repr(el.GetLength()) + " and Radius = "
                        + repr(el.GetMeanRadius()) + ".")
                    indexUniqueBranches += 1
        else:
            el = Element(Id = indexUniqueBranches)
            el.SetMeanRadius(ComputeGroupRadius(centerline, i))
            el.SetLength(ComputeGroupLength(centerline, i))
            groupIdsArray = centerline.GetCellData().GetArray(GROUPIDSARRAYNAME)
            numberOfCells = centerline.GetNumberOfCells()
            for k in range(0, numberOfCells):
                groupId = groupIdsArray.GetValue(k)
                if groupId != i:
                    continue
                x0 = [0.0, 0.0, 0.0]
                x1 = [0.0, 0.0, 0.0]
                npts = centerline.GetCell(k).GetPoints().GetNumberOfPoints()
                centerline.GetCell(k).GetPoints().GetPoint(npts - 1, x1)
                x1List.append(x1)
                centerline.GetCell(k).GetPoints().GetPoint(0, x0)
                x0List.append(x0)
                VtkCellIdList.append(k)
                VtkGroupIdList.append(i)
            uniqueX0List = [list(x) for x in set(tuple(x) for x in x0List)]
            uniqueX1List = [list(x) for x in set(tuple(x) for x in x1List)]
            el.SetInOutPointsCoordinates(uniqueX0List, uniqueX1List)
            el.SetVtkGroupIdList(VtkGroupIdList)
            el.SetVtkCellIdList(VtkCellIdList)
            network.AddElement(el)
            verboseprint("> Edge Id " + repr(indexUniqueBranches) 
                + ", Length = " + repr(el.GetLength()) + " and Radius = "
                + repr(el.GetMeanRadius()) + ".")
            indexUniqueBranches += 1
    verboseprint("> ")

    if isConnectivityNeeded:
        minLength, maxLength = ComputeGeometricTolerance(centerline)
        ComputeConnectivity(network, minLength, verboseprint)
    SetRadiusX0(centerline, network, verboseprint)
    network.SetNetworkInletRadius(
        ComputeInletAverageRadius(centerline, 0.0, verboseprint))

    return centerline

def ComputeConnectivity(network, tolerance, verboseprint):
    '''Compute the branches connectivity in the network.

    It is assumed that the first element is the inlet. In case of several
    inlets, modification to this routine should be coded. This routine search
    for each branch, its connected branches by comparing the coordinates of 
    the distal extremity of the treated branch with the proximal extremity
    of the other branches in the network. As sometimes the coordinates do not
    perfectly match, a tolerance is applied.

    TODO: use simply vtkPolyDataConnectivityFilter to do the job...?

    '''
    # A few cases were bordeline, hence the factor.
    tol = 5.0 * tolerance
    # Initialization of the first branch.
    network.elements[0].SetIfInlet(True)   
    network.elements[0].SetInOutPointsIds(1, 2)
    # Test the end node (x1) of a branch with the start node
    # of the other branches (x0).
    for treatedBranch in network.elements:
        atLeastOneFound = False
        for x1 in treatedBranch.GetOutPointsx1():
            for otherBranch in network.elements:
                if otherBranch.GetId() == treatedBranch.GetId():
                    continue
                for x0 in otherBranch.GetInPointsx0():
                    if vtk.vtkMath.Distance2BetweenPoints(x0,x1) < tol:
                        if vtk.vtkMath.Distance2BetweenPoints(x0,x1) > tolerance:
                            print('WARNING: POTENTIAL CONNECTIVITY ISSUE. ')
                            print('         A distance between connected points is suspicious.') 
                            print('         The segment(s) of CELL ID VTK ' + repr(treatedBranch.GetVtkCellIdList()))
                            print('         and, the segment(s) of CELL ID VTK ' + repr(otherBranch.GetVtkCellIdList()) + ' will be considered') 
                            print('         as connected. Please check if this action was expected.')
                            print()
                        otherBranch.SetInOutPointsIds(
                            treatedBranch.GetOutPointsx1Id(), 
                            otherBranch.GetId() + 2)
                        otherBranch.SetBehindSegment(treatedBranch.GetId())
                        treatedBranch.SetFrontSegment(otherBranch.GetId())
                        atLeastOneFound = True
        
        if not(atLeastOneFound):
            treatedBranch.SetIfOutlet(True)
            treatedBranch.SetInOutPointsIds(
                treatedBranch.GetInPointsx0Id(), 
                treatedBranch.GetId() + 2)

def SetRadiusX0(centerline, network, verboseprint):
    for branch in network.elements:
        cellID = branch.GetVtkCellIdList()
        radiusArray = centerline.GetPointData().GetArray(RADIUSARRAYNAME)
        r = radiusArray.GetComponent(centerline.GetCell(cellID[0]).GetPointId(0),0)
        branch.SetInletRadius(r)

def ComputeInletAverageRadius(centerline, desiredLength, verboseprint):
    '''Compute the inlet radius as an averaged radius.

    Computes an average radius over a certain length of the ICA. The mean
    radius is computed with the maximum inscribed sphere radius in the vessel.
    Thus, the actual radius might be underestimated if the vessel has
    an elliptical shape.

    '''
    branchId = 0
    length = 0.0
    resistance = 0.0
    radiusArray = centerline.GetPointData().GetArray(RADIUSARRAYNAME)
    npts = GetIndexCenterlineForADefinedLength(centerline, branchId, 
        desiredLength, verboseprint)
    for k in range(0, npts - 1):
        point0 = [0.0, 0.0, 0.0]
        point1 = [0.0, 0.0, 0.0]
        centerline.GetCell(branchId).GetPoints().GetPoint(k, point0)
        centerline.GetCell(branchId).GetPoints().GetPoint(k + 1, point1)
        dx = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))
        length += dx
        r = radiusArray.GetComponent(centerline.GetCell(branchId).GetPointId(k),0)
        resistance += dx / r**(4.0)
    if npts < 2:
        branchMeanRadius = radiusArray.GetComponent(
            centerline.GetCell(branchId).GetPointId(0),0)
    else:
        branchMeanRadius = (length / resistance)**(0.25)

    return branchMeanRadius

def ComputeInletAverageCrossSectionArea(centerline, desiredLength, verboseprint):
    '''Compute the inlet radius as an averaged radius.

    Computes an average radius over a certain length of the ICA. The mean
    radius is computed with the cross section area. If the number of points
    over the segment is lower than two, then the section area correspond
    to the inlet point. Note that due to a bug in VTK the very first point
    seems wrong, that is why I am taking the second point of the centerline.

    '''
    branchId = 0
    length = 0.0
    resistance = 0.0
    sectionArray = centerline.GetPointData().GetArray(SECTIONARRAYNAME)
    npts = GetIndexCenterlineForADefinedLength(centerline, branchId, 
        desiredLength, verboseprint)
    for k in range(1, npts - 1):
        point0 = [0.0, 0.0, 0.0]
        point1 = [0.0, 0.0, 0.0]
        S = sectionArray.GetComponent(centerline.GetCell(branchId).GetPointId(k),0)
        if S == 0.0:
            continue
        centerline.GetCell(branchId).GetPoints().GetPoint(k, point0)
        centerline.GetCell(branchId).GetPoints().GetPoint(k + 1, point1)
        dx = math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))
        length += dx
        r = math.sqrt(S/np.pi)
        resistance += dx / r**(4.0)
    if npts < 2:
        S = sectionArray.GetComponent(centerline.GetCell(branchId).GetPointId(1),0)
        branchMeanRadius = math.sqrt(S/np.pi)
    else:
        branchMeanRadius = (length / resistance)**(0.25)

    return branchMeanRadius

def GetIndexCenterlineForADefinedLength(centerline, branchId, desiredLength, 
    verboseprint):
    '''Get the index of the point such as the desired distance between the index and 
    the beginning of the branch is reached. '''

    length = 0.0
    done = False
    npts = centerline.GetCell(branchId).GetPoints().GetNumberOfPoints()
    for k in range(0, npts - 1):
        point0 = [0.0, 0.0, 0.0]
        point1 = [0.0, 0.0, 0.0]
        centerline.GetCell(branchId).GetPoints().GetPoint(k, point0)
        centerline.GetCell(branchId).GetPoints().GetPoint(k + 1, point1)
        length += math.sqrt(vtk.vtkMath.Distance2BetweenPoints(point0,point1))
        if length >= desiredLength:
            index = k
            done = True
            break
    if not(done):
        try:
            index = k
        except:
            index = None

    return index

def GetListProbePoints(centerline, network, verboseprint):
    '''Get points on a centerline spaced by the inlet diameter. '''

    diameterInlet = 2.0*network.GetNetworkInletRadius()
    pointsList = []
    for el in network.elements:
        if el.IsBlanked():
            continue
        branchId = el.GetVtkCellIdList()[0]
        npts = centerline.GetCell(branchId).GetPoints().GetNumberOfPoints()
        done = False
        ctr = 0
        while not(done):
            x = [0.0, 0.0, 0.0]
            desiredLength = float(ctr)*diameterInlet
            idx = GetIndexCenterlineForADefinedLength(centerline, branchId, desiredLength, 
                    verboseprint)
            if idx is None:
                break
            centerline.GetCell(branchId).GetPoints().GetPoint(idx, x)
            pointsList.append(x) 
            ctr += 1
            if idx == npts - 2:
                done = True
                
    return(pointsList)


class Element(object):

    def __init__(self, Id):
        self.Id = Id
        self.vtkGroupIdList = []
        self.vtkCellIdList = []
        self.length = 0.0
        self.inletRadius = 0.0
        self.meanInletRadius = 0.0
        self.outletRadius = 0.0
        self.meanRadius = 0.0
        self.blanking = 0
        self.inlet = False
        self.outlet = False
        self.x0 = []
        self.x1 = []
        self.x0Id = -1
        self.x1Id = -1
        self.behindSegment = None
        self.frontSegment = None
        self.alpha = 1.0
        self.beta = 0.0
        self.gamma = 0.0

    def SetLength(self, length):
        self.length = length

    def SetAlpha(self, alpha):
        self.alpha = alpha

    def SetBeta(self, beta):
        self.beta = beta

    def SetGamma(self, gamma):
        self.gamma = gamma

    def SetBehindSegment(self, id):
        self.behindSegment = id

    def SetFrontSegment(self, id):
        self.frontSegment = id

    def SetVtkGroupIdList(self, VtkGroupIdList):
        self.vtkGroupIdList = VtkGroupIdList

    def SetVtkCellIdList(self, vtkCellIdList):
        self.vtkCellIdList = vtkCellIdList

    def SetMeanRadius(self, radius):
        self.meanRadius = radius

    def SetInletRadius(self, radius):
        self.inletRadius = radius

    def SetOutletRadius(self, radius):
        self.outletRadius = radius

    def SetBlanking(self, blanking):
        self.blanking = blanking

    def SetIfInlet(self, inlet):
        self.inlet = inlet

    def SetIfOutlet(self, outlet):
        self.outlet = outlet

    def SetInOutPointsCoordinates(self, x0, x1):
        self.x0 = x0
        self.x1 = x1

    def SetInOutPointsIds(self, x0Id, x1Id):
        self.x0Id = x0Id
        self.x1Id = x1Id

    def GetId(self):
        return self.Id

    def GetBehindSegment(self):
        return self.behindSegment

    def GetFrontSegment(self):
        return self.frontSegment

    def GetVtkCellIdList(self):
        return self.vtkCellIdList

    def GetVtkGroupIdList(self):
        return self.vtkGroupIdList

    def GetLength(self):
        return self.length

    def GetAlpha(self):
        return self.alpha

    def GetBeta(self):
        return self.beta

    def GetMeanRadius(self):
        return self.meanRadius

    def GetInletRadius(self):
        return self.inletRadius

    def GetMeanArea(self):
        meanArea = np.pi*(self.meanRadius**2.0)
        return meanArea

    def GetInPointsx0(self):
        return self.x0

    def GetOutPointsx1(self):
        return self.x1

    def GetInPointsx0Id(self):
        return self.x0Id

    def GetOutPointsx1Id(self):
        return self.x1Id

    def IsBlanked(self):
        '''Returns True if the element is part of a branch division.'''
        return self.blanking

    def IsAnInlet(self):
        return self.inlet

    def IsAnOutlet(self):
        return self.outlet


class Network(object):

    def __init__(self):
        self.elements = []
        self.numberOfElements = 0
        self.numberOfBifurcations = 0
        self.numberOfOutlets = 0
        self.networkInletRadius = 0.0

    def AddElement(self, x):
        self.elements.append(x)
        self.numberOfElements += 1

    def SetNetworkInletRadius(self, radius):
        self.networkInletRadius = radius
    
    def GetNumberOfElements(self):
        return self.numberOfElements

    def GetNumberOfBifBranches(self):
        nBlanked = 0
        for el in self.elements:
            nBlanked += el.IsBlanked()
        return nBlanked

    def GetNumberOfInlet(self):
        nInlet = 0
        for el in self.elements:
            nInlet += el.IsAnInlet()
        return nInlet

    def GetNumberOfOutlet(self):
        nOutlet = 0
        for el in self.elements:
            nOutlet += el.IsAnOutlet()
        return nOutlet

    def GetNetworkInletRadius(self):
        return self.networkInletRadius

        
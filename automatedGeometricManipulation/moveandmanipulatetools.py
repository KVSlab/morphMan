from common import *
from argparse import ArgumentParser
from os import path, listdir
from subprocess import STDOUT, check_output
from IPython import embed

import matlab.engine
import operator
import sys
import math
import numpy.linalg as la
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema, resample
from scipy.ndimage.filters import gaussian_filter as gauss
from mpl_toolkits.mplot3d import Axes3D

# Local import
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from matplotlib import rc, rcParams
from vmtk import vtkvmtk, vmtkscripts

def vmtkSmoother(surface, method, iterations=600):
    smoother= vmtkscripts.vmtkSurfaceSmoothing()
    smoother.Surface = surface
    smoother.NumberOfIterations = iterations
    smoother.Method = method
    smoother.Execute()
    surface = smoother.Surface

    return surface

def get_points(data, key, bif=False):
    if bif:
        div_points = np.asarray([data[0][key], data[1][key]])
    else:
        div_points = np.asarray([data["bif"][key], data[0][key], data[1][key]])

    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)

    return points, div_points
def get_clipping_points(dirpath):
    particles = path.join(dirpath, "optimal_bog.particles")
    all_points = np.loadtxt(particles)
    clipping_points = all_points[2:] 
    return clipping_points

def sort_centerlines(centerlines_complete):
    # Find longest line:
    lines = []
    n = centerlines_complete.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines_complete, i))

    longest = [lines[0]]
    lenlong = get_curvilinear_coordinate(longest[0])
    for i in range(1,n):
        tmplong = get_curvilinear_coordinate(lines[i])
        if len(tmplong) > len(lenlong):
            lenlong = tmplong
            longest.insert(0, lines[i])
        else: 
            longest.append(lines[i])

    centerlines_in_order = merge_data(longest)
    return centerlines_in_order


def connect_line(line):
    
    # Create edges between new_centerline points
    pts = vtk.vtkPoints()
    for i in range((line.GetNumberOfCells())):
        pts.InsertNextPoint(line.GetPoint(i))

    lines = vtk.vtkCellArray()
    for i in range((line.GetNumberOfCells()) - 1):
        newline = vtk.vtkLine()
        newline.GetPointIds().SetId(0, i)
        newline.GetPointIds().SetId(1, i + 1)
        lines.InsertNextCell(newline)

    line = vtk.vtkPolyData()
    line.SetPoints(pts)
    line.SetLines(lines)
    return line


def move_line_horo(patch_cl, ID1, ID2, dx_p1, clip=False, eye=False, side=None):

    centerline_loc = get_locator(patch_cl)
    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if clip == True:
        if eye == True:
            l1 = ExtractSingleLine(patch_cl, 0)
            l2 = ExtractSingleLine(patch_cl, 1)
            l3 = ExtractSingleLine(patch_cl, 2)
            test_cl = merge_data([l1,l2,l3])

            ID1 = 0
            ID2 = len(get_curvilinear_coordinate(test_cl))
            idmid = int( (ID1 + ID2)/2.)
        else:
            ID1 = 0
            ID2 = len(get_curvilinear_coordinate(patch_cl))
            idmid = int( (ID1 + ID2)/2.)
        
        for p in range(patch_cl.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(patch_cl.GetPoint(p))

            if cl_id < idmid:
                dist = dx_p1 * (idmid**2-cl_id**2) / (idmid**2-ID1**2)
            else:
                if cl_id <= (ID2 - 1):
                    dist = -dx_p1 * (cl_id-idmid)**(0.5) / (ID2-idmid)**(0.5)
                else: 
                    locator = get_locator(test_cl)
                    pp = patch_cl.GetPoint(cl_id)
                    id_main = locator.FindClosestPoint(pp)
                    dist = -dx_p1 * (id_main-idmid)**(0.5) / (id_main-idmid)**(0.5)

            patch_point = np.asarray(patch_cl.GetPoint(p))

            if la.norm(patch_point) > 0.1:
                points.InsertNextPoint(patch_point + dist)
                verts.InsertNextCell(1)
                verts.InsertCellPoint(p)

    else:
        for p in range(patch_cl.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(patch_cl.GetPoint(p))

            if side == "right":
                dist = -dx_p1
            elif side == "left":
                dist = dx_p1


            points.InsertNextPoint(np.asarray(patch_cl.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(patch_cl.GetPointData().GetArray(radiusArrayName))

    return newline



def spline_line(line, get_curv=False, isline=False, nknots=50):
    
    if isline == False:
        # Create edges between new_centerline points
        pts = vtk.vtkPoints()
        for i in range((line.GetNumberOfCells())):
            pts.InsertNextPoint(line.GetPoint(i))

        lines = vtk.vtkCellArray()
        for i in range((line.GetNumberOfCells()) - 1):
            newline = vtk.vtkLine()
            newline.GetPointIds().SetId(0, i)
            newline.GetPointIds().SetId(1, i + 1)
            lines.InsertNextCell(newline)

        line = vtk.vtkPolyData()
        line.SetPoints(pts)
        line.SetLines(lines)

    # Collect data from centerline
    data = np.zeros((line.GetNumberOfPoints(), 3))

    curv_coor = get_curvilinear_coordinate(line)
    for i in range(data.shape[0]):
        data[i,:] = line.GetPoint(i)

    t = np.linspace(curv_coor[0], curv_coor[-1], nknots+2)[1:-1]
    fx = splrep(curv_coor, data[:,0], k=4, t=t)
    fy = splrep(curv_coor, data[:,1], k=4, t=t)
    fz = splrep(curv_coor, data[:,2], k=4, t=t)

    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    data = np.zeros((len(curv_coor), 3))
    data[:,0] = fx_
    data[:,1] = fy_
    data[:,2] = fz_

    header = ["X", "Y", "Z"]
    line = data_to_vtkPolyData(data, header)
   
    # Let vmtk compute curve attributes
    line = CenterlineAttribiutes(line, smooth=False)

    if get_curv == True:
        # Analytical curvature
        dlsfx = splev(curv_coor, fx, der=1)
        dlsfy = splev(curv_coor, fy, der=1)
        dlsfz = splev(curv_coor, fz, der=1)

        ddlsfx = splev(curv_coor, fx, der=2)
        ddlsfy = splev(curv_coor, fy, der=2)
        ddlsfz = splev(curv_coor, fz, der=2)

        C1xC2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
        C1xC2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
        C1xC2_3 = ddlsfy * dlsfx - ddlsfx * dlsfy

        curvature = np.sqrt(C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2) / \
                            (dlsfx**2 + dlsfy**2 + dlsfz**2)**1.5

        return line, curvature
    else:
        return line



def move_dots_vertical(line, dx, ID1_0, clip_ID=None, eye=False):

    centerline_loc = get_locator(line)

    newline = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    if eye == True:
        l1 = ExtractSingleLine(line, 0)
        l2 = ExtractSingleLine(line, 1)
        l3 = ExtractSingleLine(line, 2)
        test_cl = merge_data([l1,l2,l3])

        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(test_cl))
        IDmid = int( (ID1 + ID2)/2.)
        I1 = I2 = IDmid


        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))
            if cl_id < I1:
                I1 = cl_id

            if cl_id > I2 and cl_id <= (ID2-4):
                I2 = cl_id

        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))

            if cl_id <= (ID2-4):
                dist = 4 * dx * (cl_id - I1)*(I2 - cl_id) / ( I2 - I1)**2
            else:
                cl_id = clip_ID - ID1_0 + int((ID2 - (clip_ID - ID1_0) )*0.4)
                dist = 4 * dx * (cl_id - ID1)*(ID2 - cl_id) / ( ID2 - ID1)**2

            points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    else:
        ID1 = 0
        ID2 = len(get_curvilinear_coordinate(line))
    
        for p in range(line.GetNumberOfPoints()):
            cl_id = centerline_loc.FindClosestPoint(line.GetPoint(p))

            dist = 4 * dx * (cl_id - ID1)*(ID2 - cl_id) / ( ID2 - ID1)**2

            points.InsertNextPoint(np.asarray(line.GetPoint(p)) + dist)
            verts.InsertNextCell(1)
            verts.InsertCellPoint(p)

    newline.SetPoints(points)
    newline.SetVerts(verts)
    newline.GetPointData().AddArray(line.GetPointData().GetArray(radiusArrayName))
    return newline


def movePerp(xx,yy,zz,n, P, Z, alpha):
    p1 = P[0]
    p2 = P[1]

    # Find midpoint and point furthest away
    dist = []
    for z in Z:
        d = la.norm( np.cross((z - p1),(z - p2))) / la.norm(p2 - p1)
        dist.append(d)

    D_id = dist.index(max(dist))
    D = max(dist)
    z_m = Z[D_id]

    # Vector from line to Z_max and projection onto plane
    v = (z_m - p2) - (z_m - p2).dot(p2-p1)*(p2-p1) / la.norm(p2-p1)**2
    PV = v - v.dot(n)*n

    # Find distances 
    P1 =  (z_m - p1).dot(p2-p1)*(p2-p1) / la.norm(p2-p1)**2
    P1 = p1 + P1
    V = P1 + v
    PV1 = P1 + PV

    # Move points
    dZ = []
    for i in range(len(Z)):
        dz = np.array(PV) * dist[i] / D * alpha
        dZ.append( Z[i] + dz)

    return dZ, (PV1 - P1)*alpha

def movePara(xx,yy,zz,n, P,Z, beta, case ):
    p1 = P[0]
    p2 = P[1]
    
    # Find normal from q
    q = [(p1 + p2)/2.]
    qp2 = np.array(p2 - q)[0]
    qp1 = np.array(p1 - q)[0]
    s = q[0] - np.cross(qp2, n)*3

    # Split points based on orientation
    # to q normal
    Z_p = []
    Z_p_dist = []
    Z_m = []
    Z_m_dist = []
    for z in Z:
        d = la.norm( np.cross((z - s),(z - q[0]))) / la.norm(q[0] - s)
        c = np.cross(s-q, z-q)
        if c[0][0] >= 0:
            Z_p.append(z)
            Z_p_dist.append(d)
        else:
            Z_m.append(z)
            Z_m_dist.append(d)
    
    # Move points
    D = la.norm(p1-p2) / 2.
    dZ = []
    dZ.append(p1 + qp1*beta)
    for i in range(len(Z_p)):
        dz = qp1 * Z_p_dist[i] / D * beta
        dZ.append(Z_p[i] + dz)

    for i in range(len(Z_m)):
        dz = qp2 * Z_m_dist[i] / D * beta
        dZ.append(Z_m[i] + dz)
    dZ.append(p2 + qp2*beta)
    
    # Check if moved in right direction
    d_0 = la.norm( np.cross(Z_p[0] - q, Z_p[0] - s) ) / la.norm(s-q)
    d_1 = la.norm( np.cross(dZ[1] - q, dZ[1] - s) ) / la.norm(s-q)
    if d_1 < d_0:
        # Split points based on orientation
        # to q normal
        Z_p = []
        Z_p_dist = []
        Z_m = []
        Z_m_dist = []
        for z in Z:
            d = la.norm( np.cross((z - s),(z - q[0]))) / la.norm(q[0] - s)
            c = -np.cross(s-q, z-q)
            if c[0][0] >= 0:
                Z_p.append(z)
                Z_p_dist.append(d)
            else:
                Z_m.append(z)
                Z_m_dist.append(d)

        # Move points
        D = la.norm(p1-p2) / 2.
        dZ = []
        dZ.append(p1 + qp1*beta)
        for i in range(len(Z_p)):
            dz = qp1 * Z_p_dist[i] / D * beta
            dZ.append(Z_p[i] + dz)

        for i in range(len(Z_m)):
            dz = qp2 * Z_m_dist[i] / D * beta
            dZ.append(Z_m[i] + dz)
        dZ.append(p2 + qp2*beta)


    zpid = Z_p_dist.index(min(Z_p_dist))
    zp_min = Z_p[zpid]
    zmid = Z_m_dist.index(min(Z_m_dist))
    zm_min = Z_m[zmid]
    return dZ, zp_min, zm_min

def get_clipped_curve(centerlines_in_order, clipping_points, end=False): 
    # Extract clipped segment
    start_p = vtk.vtkPoints()
    end_p   = vtk.vtkPoints()
    line    = ExtractSingleLine(centerlines_in_order, 0)
    start   = line.GetPoint(0)
    endId   = len(get_curvilinear_coordinate(line)) -2
    end     = line.GetPoint(endId)
    for p1, p2 in zip([start, clipping_points[0]], [clipping_points[1], end]):
        start_p.InsertNextPoint(p1)
        end_p.InsertNextPoint(p2)

    clip_end = CreateParentArteryPatches(line, end_p,True)

    if end == True:
        cut_end = vtk.vtkPoints()
        cut_end.InsertNextPoint(start)
        cut_end.InsertNextPoint(clipping_points[1])
        return CreateParentArteryPatches(centerlines_in_order, cut_end, True)
    else: 
        clipped_curve = CreateParentArteryPatches(clip_end, start_p,True)

        return clipped_curve



def bestPlane(Z, P):
    # Defined matrices
    A = Z
    C = P
    b = np.ones(len(A))
    d = np.ones(len(C))

    # Create complete matrix
    ATA = np.transpose(A).dot(A)
    M0 = np.c_[ATA,np.transpose(C)]
    M1 = np.c_[C, np.zeros( (len(C), len(C)))]
    M = np.r_[M0, M1]
    Y = np.r_[np.transpose(A).dot(b),d]
    
    # Solve system
    x = la.solve(M,Y)
    a = x[0]
    b = x[1]
    c = x[2]
    n = np.array([a,b,c])
    n = n / la.norm(n)

    # Define plane
    xmin = min(Z, key=operator.itemgetter(1))[0] - 4 
    xmax = max(Z, key=operator.itemgetter(1))[0] + 4
    ymin = min(Z, key=operator.itemgetter(1))[1] - 4 
    ymax = max(Z, key=operator.itemgetter(1))[1] + 4
    xx,yy = np.meshgrid(np.linspace(xmin,xmax,15),np.linspace(ymin,ymax,15))
    zz = (1 - a*xx - b*yy ) / float(c)
    return xx,yy,zz,n


def get_spline_points(line, param,growDir,case, clip_points, voronoi):
    """ 
    Pick n uniformly selected points along the 
    centerline from point P1 to P2
    And move them
    """
    locator = get_locator(line)
    p1      = clip_points.GetPoint(0)
    p2      = clip_points.GetPoint(1)
    ID1     = locator.FindClosestPoint(p1)
    ID2     = locator.FindClosestPoint(p2)
    ID_mid = int((ID1 + ID2) / 2.)
    P = [p1,p2]
    
    # Select n uniformly spaced points 
    n = 10
    points = []
    ids = np.zeros(n)
    dx = 1 / (n + 1.)
    for i in range(1,n+1):
        ID = int(ID1 + (ID2 - ID1) * i * dx)
        ids[i-1] = ID
        p = line.GetPoints().GetPoint(ID)
        points.append(np.array([p[0], p[1], p[2]]))

    for i in range(len(P)):
        P[i] = np.array([P[i][0], P[i][1], P[i][2]])

    xx,yy,zz,n = bestPlane(points,P)

    if growDir == "vertical":
        dz, dx = movePerp(xx,yy,zz,n,P, points, param)
        return dz, ids, dx
    else: 
        dz,zp,zm = movePara(xx,yy,zz,n,P, points, param, case)
        return dz, ids



def find_ophthalmic_artery(centerlines, clip_pts):

    # Extract lines:
    lines = []
    n = centerlines.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines, i))

    longest = lines[0]
    tol = 0.40
    
    # Find start and stop IDs along clipped curve
    locator_longest = get_locator(longest)
    p1      = clip_pts[0]
    p2      = clip_pts[1]
    ID1     = locator_longest.FindClosestPoint(p1)
    ID2     = locator_longest.FindClosestPoint(p2)
    
    eye = False
    index = 1
    for line in lines[1:]:
        locator_checkline = get_locator(line)
        len_check = len(get_curvilinear_coordinate(line)) 
        if len_check < ID2:
            IDStop = len_check - 1
        else:
            IDStop = ID2

        for i in np.arange(ID1, IDStop):
            p_eye   = np.asarray(line.GetPoint(i))
            p_cl    = np.asarray(longest.GetPoint(i))
            dist = la.norm(p_eye - p_cl)
            if  dist > tol:
                clip_ID = i
                eye = True
                eyeline = lines[index]
                del lines[index]
                centerlines = merge_data(lines)
                break
        if eye:
            break
        index += 1
    
    if eye:
        return eye, clip_ID, centerlines, eyeline
    else:
        return eye, None, centerlines, None

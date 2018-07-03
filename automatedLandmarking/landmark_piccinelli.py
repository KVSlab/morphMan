#!/usr/bin/env python
import matplotlib.pyplot as plt
import matlab.engine
import numpy as np
import numpy.linalg as la
import sys
import math

from scipy import signal
from common import *
from os import path, listdir
from scipy.interpolate import splrep, splev
from scipy.signal import argrelextrema 
from scipy.ndimage.filters import gaussian_filter
from subprocess import check_output, STDOUT
from IPython import embed

def landmarking_piccinelli(centerline, folder, step=None, it=None, fac_c=None, fac_t=None, nknots=None):
    """ Takes the file path to a centerline for further resampling followed by smoothing
    using vmtkcenterlinegeometry. Landmarking is based on Piccinelli et al. (2011)
    """
    print "Case: %s" % folder

    if step != None:
        centerline = CenterlineResampling(centerline, step)

    if nknots != None:
        line, max_point_ids, min_point_ids = splineCenterline(centerline,nknots)
        curvature = get_array("Curvature", line)
        torsion = get_array("Torsion", line)
        torsion_smooth = gaussian_filter(torsion , 2)

        min_tor = list(argrelextrema(torsion_smooth, np.less)[0])
        max_tor = list(argrelextrema(torsion_smooth, np.greater)[0])
        extreme_tor = sorted(min_tor + max_tor)

    else:
        line = CenterlineGeometry(centerline, outputsmoothed=0, factor=fac_c, iterations=it)
        line_tor = CenterlineGeometry(centerline, outputsmoothed=0, factor=fac_t, iterations=it)

        curvature = get_array("Curvature", line)
        torsion = get_array("Torsion", line_tor)
        torsion_smooth = gaussian_filter(torsion , 2)

        max_point_ids = list(argrelextrema(curvature, np.greater)[0])
        min_tor = list(argrelextrema(torsion_smooth, np.less)[0])
        max_tor = list(argrelextrema(torsion_smooth, np.greater)[0])
        extreme_tor = sorted(min_tor + max_tor)

    # Extract local curvature minimums
    length = get_curvilinear_coordinate(line)
    
    # Remove points too close to the end of the siphon
    for i in max_point_ids:
        if length[i] in length[-5:]:
            max_point_ids.remove(i)
    
    # Remove curvature peaks too close to each other
    threshold = 30
    pointdistance= []
    for i in range(0, len(max_point_ids)-1):
        pointdistance.append(max_point_ids[i+1] - max_point_ids[i])

    for i, d in enumerate(pointdistance):
        if pointdistance[i] < threshold:
            curv1 = curvature[i]
            curv2 = curvature[i+1]
            if curv1 > curv2:

                max_point_ids[i+1] = None
            else:
                max_point_ids[i] = None

    max_point_ids = [item for item in max_point_ids if item != None]

    # Remove torsion saddle points 
    tor_extreme = extreme_tor[:]
    for i, d in enumerate(extreme_tor):
        tor1 = torsion_smooth[d][0]
        tol = 0.3
        if abs(tor1) < tol:
            tor_extreme.remove(d)

    # Define bend interfaces based on Piccinelli et al.
    def find_interface():
        interface = {}
        k=0
        for c in max_point_ids:
            for i in range(0, len(tor_extreme) - 1):
                if tor_extreme[i] < c and c < tor_extreme[i+1]:
                    interface["bend%s" % (k+1)] = np.array([tor_extreme[i]])
                    k+=1
                    interface["bend%s" % (k+1)] = np.array([tor_extreme[i+1]])
                    k+=1
        
        return interface
    
    # Compute and extract interface points
    interfaces = find_interface()
    landmarks = {}
    for k, v in interfaces.iteritems():
        landmarks[k] = line.GetPoints().GetPoint(int(v))

    if landmarks is not None:
        writeParameters(landmarks, folder)
        createParticles(folder)

def splineCenterline(line, nknots, step=None):
    # Allocate data structure to store centerline points
    data = np.zeros((line.GetNumberOfPoints(), 3))

    # Collect data from centerline
    curv_coor = get_curvilinear_coordinate(line)
    for i in range(data.shape[0]):
        data[i,:] = line.GetPoints().GetPoint(i)

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
    length = get_curvilinear_coordinate(line)
    
    # Compute curvature from the 'exact' spline to get a robust way of
    # finding max / min points on the centerline

    dlsfx = splev(curv_coor, fx, der=1)
    dlsfy = splev(curv_coor, fy, der=1)
    dlsfz = splev(curv_coor, fz, der=1)

    ddlsfx = splev(curv_coor, fx, der=2)
    ddlsfy = splev(curv_coor, fy, der=2)
    ddlsfz = splev(curv_coor, fz, der=2)

    C1xC2_1 = ddlsfz * dlsfy - ddlsfy * dlsfz
    C1xC2_2 = ddlsfx * dlsfz - ddlsfz * dlsfx
    C1xC2_3 = ddlsfy * dlsfx - ddlsfx * dlsfy

    curvature_ = np.sqrt(C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2) / \
                        (dlsfx**2 + dlsfy**2 + dlsfz**2)**1.5

    max_point_ids = list(argrelextrema(curvature_, np.greater)[0])
    min_point_ids = list(argrelextrema(curvature_, np.less)[0])
    
    # TODO: Replace the locator with curv_coor = length
    locator = get_locator(line)

    min_points = [[fx_[i], fy_[i], fz_[i]] for i in min_point_ids]
    max_points = [[fx_[i], fy_[i], fz_[i]] for i in max_point_ids]
    min_point_ids = []
    max_point_ids = []

    for point_min, point_max in zip(min_points, max_points):
        min_point_ids.append(locator.FindClosestPoint(point_min))
        max_point_ids.append(locator.FindClosestPoint(point_max))

    # The ParallelTransportNormals and the FrenetTangent is not orthonormal
    # (but close) from vmtk. Using the GramSchmidt proses gives E2 and fixes
    # the non-orthogonality
    E1 = get_array("ParallelTransportNormals", line, k=3)
    T = get_array("FrenetTangent", line, k=3)
    E2 = np.zeros((E1.shape[0], 3))

    for i in range(E1.shape[0]):
        V = np.eye(3)
        V[:, 0] = T[i,:]
        V[:, 1] = E1[i,:]
        V = GramSchmidt(V)

        E1[i,:] = V[:,1]
        E2[i,:] = V[:,2]

    # Compute k_1, k_2 furfilling T' = curv(s)*N(s) = k_1(s)*E_1(s) + k_2(s)*E_2(s).
    # This is simply a change of basis for the curvature vector N. The room
    # of k_1 and k_2 can be used to express both the curvature and the
    # torsion.
    curvature = get_array("Curvature", line)
    # Compute a improved torsion that does not include all the noice from
    # vmtk. Here the analytical representation of the spline is accessebole.
    length = get_curvilinear_coordinate(line)
    dddlsfx = splev(length, fx, der=3)
    dddlsfy = splev(length, fy, der=3)
    dddlsfz = splev(length, fz, der=3)

    torsion_spline = (dddlsfx*C1xC2_1 + dddlsfy*C1xC2_2 + dddlsfz * C1xC2_3) / \
                        (C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2)
    torsion_array = create_vtk_array(torsion_spline, "Torsion")
    line.GetPointData().AddArray(torsion_array)
    
    curvature_ = np.sqrt(C1xC2_1**2 + C1xC2_2**2 + C1xC2_3**2) / \
                        (dlsfx**2 + dlsfy**2 + dlsfz**2)**1.5
    curvature_[0] = curvature[0]
    curvature_[-1] = curvature[-1]

    curvature_array = create_vtk_array(curvature_, "Curvature")
    line.GetPointData().AddArray(curvature_array)

    return line, max_point_ids, min_point_ids


def createParticles(folder):
    # Create a file with points where bends are located and 
    # remove points from manifest
    manifest = folder + "/manifest.txt"
    filename = folder + "/piccinelli.particles" 

    output = open(filename, "w")
    mani= open(manifest, "r")
    lines = mani.readlines()
    mani.close()
    mani = open(manifest, "w")

    for line in lines:
        if line.split()[1][0] == "(":
            coord = line.split()[1:]
            point = "%s %s %s" % (coord[0][1:-1], coord[1][:-1], coord[2][:-1])
            output.write(point + "\n")
        else:
            mani.write(line)

    output.close()
    mani.close()


def main(folder, method):
    case = folder.split("/")[-1]
    
    line = ExtractCarotidSiphon(folder)

    if method == "vmtk":
        # Input parameters
        # VMTK:
        resamp_step = 0.1
        it = 100
        fac_curvature = 1.5
        fac_torsion = 1.2
        landmarking_piccinelli(line, folder, step=resamp_step, it=it, fac_c=fac_curvature, fac_t=fac_torsion)
    elif method == "splines":
        # B - Splines:
        resamp_step = None
        nknots = get_knots(line)
        landmarking_piccinelli(line, folder, nknots=nknots, step=resamp_step)
    

def get_knots(centerline):
    L = centerline.GetLength()
    return np.ceil( L / 20. + 26. / 3. )


if __name__ == '__main__':
    maindir = '/home/henrik/master_uio/cases/'
    method = "vmtk" # ["vmtk", "splines"]

    for folder in sorted(listdir(maindir)):
        dirpath = path.join(maindir, folder)
        if path.isdir(dirpath):
            main(dirpath, method)


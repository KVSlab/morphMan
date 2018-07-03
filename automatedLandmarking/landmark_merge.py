# -*- coding: utf-8 -*-
#!/usr/bin/env python
import matplotlib.pyplot as plt
import matlab.engine
import numpy as np
import sys
import math
import numpy.linalg as la

from common import *
from patchandinterpolatecenterlines import *
from os import path, listdir
from scipy.interpolate import splrep, splev
from scipy.signal import argrelextrema, gaussian, resample
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage.filters import gaussian_filter as gauss
from subprocess import check_output, STDOUT
from IPython import embed


def writeFile(knots, x, lenth):
    splinefold = "/home/hkjeldsberg/master_uio/b_spline" 
    knotfile = splinefold + "/knots.txt"
    xfile = splinefold + "/xvalues.txt" 
    lenfile = splinefold + "/lenvalues.txt"

    knotout = open(knotfile, "w")
    lenout = open(lenfile, "w")
    xout = open(xfile, "w")

    for k in knots:
        knotout.write("%s\n" % k)

    for x_ in x:
        xout.write("%s\n" % x_)

    for l_ in lenth:
        lenout.write("%s\n" % l_)

    xout.close()
    knotout.close()
    lenout.close()



def splineCenterline(line, smooth, nknots, step=None):
    # Allocate data structure to store centerline points
    #print "Resampling line.."
    #line = CenterlineResampling(line, length=step)
    nknots = 20
    data = np.zeros((line.GetNumberOfPoints(), 3))

    # Collect data from centerline
    for i in range(data.shape[0]):
        curv_coor = get_curvilinear_coordinate(line)

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
    
    #xxx = data[:,0].copy()
    #writeFile(t, xxx, length)
    #sys.exit(1)
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
    N = get_array("FrenetNormal", line, k=3)
    curvature = get_array("Curvature", line)
    plt.plot(curvature, label="VMTK")

    k2 = (curvature.T * (E1[:,1]*N[:,0] - N[:,1]*E1[:,0]) / \
                        (E2[:,1]*E1[:,0] - E2[:,0]*E1[:,1]))[0]
    k1 = (-(curvature.T * N[:,0] + k2*E2[:,0]) / E1[:,0])[0]

    for k in [(k1, "k1"), (k2, "k2")]:
        k_array = create_vtk_array(k[0], k[1])
        line.GetPointData().AddArray(k_array)

    # Compute a improved torsion that does not include all the noice from
    # vmtk. Here the analytical representation of the spline is accessebole.
    temp_tor=get_array("Torsion", line)
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


# Resampling of centerline 

def CenterlineResampling(line, length, filename=None):

    if filename is None:
        filename = "tmp_re.vtp"
        WritePolyData(line, filename)

    command = ('vmtk vmtkcenterlineresampling -ifile %s -ofile %s -length %s') % (filename, filename, length)

    a = check_output(command, stderr=STDOUT, shell=True)

    status, text = success(a)

    if not status:
        print "Something went wrong when resampling the centerline"
        print text
        sys.exit(1)

    line = ReadPolyData(filename)
    
    return line 


def CenterlineGeometry(line, outputsmoothed, factor, iterations, filename=None):

    if filename is None:
        filename = "tmp_geo.vtp"
        WritePolyData(line, filename)

    command = ('vmtk vmtkcenterlinegeometry -ifile %s -ofile %s -smoothing 1 -outputsmoothed %s -iterations %s -factor %s -curvaturearray %s -torsionarray %s') % (filename, filename, outputsmoothed, iterations, factor, "Curvature", "Torsion")

    a = check_output(command, stderr=STDOUT, shell=True)

    status, text = success(a)

    if not status:
        print "Something went wrong when computing the geometry of the centerline "
        print text
        sys.exit(1)

    line = ReadPolyData(filename)
    
    return line 

def vizualize(length, curvature, curv_len, tor_len, min_curv, max_torsion, X, Y, Z, max_point_ids,max_coronal_ids):
        from matplotlib.pyplot import plot, hold, show, xlabel, ylabel, legend, axis
        
        # Plot curvature
        plot(length, curvature)
        plot(curv_len, min_curv, "ro")
        axis([0,length[-1], 0, 0.6])
        hold("on")
        #show()

        # Plot torsion
        plot(length, torsion)
        plot(tor_len, max_torsion,"ro")
        axis([0,length[-1], -1.5, 2.0])
        hold("on")
        #show()
        
        # Plot coordinate magnitudes
        plot(length, X, label="X")
        hold("on")
        plot(length, Y, label="Y")
        plot(length, Z, label="Z")
        legend()

        # 
        plot(length[max_coronal_ids], Z[max_coronal_ids], "g^")
        plot(length[max_coronal_ids], Z[max_coronal_ids], "g^")
        show()

        m = max_coronal_ids
        # Plot centerline
        viz(line, [[X[m], Y[m], Z[m]]])
        plot(length, curvature, [length[m], length[m]], [0,1])
        hold("on")
        plot(length[max_point_ids], curvature[max_point_ids], 'g^')
        show()

def landmarking_piccinelli(centerline, folder, showplot, resamp_step, smooth=False, nknots=25):
    """ Takes the file path to a centerline for further resampling followed by smoothing
    using vmtkcenterlinegeometry. Landmarking is based on Piccinelli et al. (2011)
    """
    print "Case: %s" % folder

    if smooth: 
        #line2, max_point_ids, min_point_ids = splineCenterline(centerline, smooth, nknots)
        line = CenterlineResampling(centerline, resamp_step)
        line = CenterlineGeometry(line, outputsmoothed=0, factor=1.0, iterations=100)
        #WritePolyData(line,"smoothline.vtp")

        curvature = get_array("Curvature", line)
        #plt.plot(curvature, label="Nofilter Curv")
        #torsion = get_array("Torsion", line)
        #plt.plot(torsion, label="No-Filter tor")
        #torsion = gaussian_filter(torsion , 5)


        #max_point_ids = list(argrelextrema(torsion, np.greater)[0])

        curvature = gaussian_filter(curvature, 5)
        min_point_ids = list(argrelextrema(curvature, np.less)[0])
        #plt.plot(torsion, label="Filter tor")
        #plt.plot(curvature, label="Curv")
        #plt.legend()
        #plt.show()

    else: 
        line, max_point_ids, min_point_ids = splineCenterline(centerline, smooth, nknots)
        torsion = get_array("Torsion", line)
        curvature = get_array("Curvature", line)

    # Find max coronal coordinate
    #value_index = Z[argrelextrema(Z, np.less)[0]].min()
    #max_coronal_ids = np.array(Z.tolist().index(value_index))
    #if abs(length[max_coronal_ids]-length[-1]) > 30:
    #    print "Sanity check failed"
    #    #sys.exit(0)
    #    return None


    # Extract local curvature minimums
    length = get_curvilinear_coordinate(line)
    
     #Remove points too close to the end of the siphon
    for i in min_point_ids:
    #for i in max_point_ids:
        if length[i] in length[-5:]:
            min_point_ids.remove(i)
            #max_point_ids.remove(i)

    curv_len = [length[i] for i in min_point_ids]
    #curv_len = [length[i] for i in max_point_ids]

    def find_interface(curv_len):
        interface = {}
        for i in range(len(curv_len)):
            min_point = min_point_ids[i]
            #min_point = max_point_ids[i]
            interface["bend%s" % (i+1)] = np.array([min_point])
        
        return interface
    
    #min_point_ids = np.array(min_point_ids)
    interfaces = find_interface(curv_len)

    landmarks = {}

    for k, v in interfaces.iteritems():
        landmarks[k] = line.GetPoints().GetPoint(int(v))

    #return landmarks
    method = "picinelli"
    smoothing = "smooth" if smooth else "spline"
    nknots = resamp_step if smooth else nknots
    if landmarks is not None:
        WritePolyData(line, path.join(folder, "%s_%s_%s.vtp" % (method, nknots, smoothing)))
        writeParameters(landmarks, folder)
        createParticles(folder, method, nknots, smoothing, resamp_step)



def spline_matlab(path, filename, init_knots, order, plot=True):
    """
    Perform Knot-free regresion spline on input centerline
    extracted from input filename
    """

    print("Computing knot free regression spline")
    #init_array_matlab = matlab.double( [i for i in init_array] )
    mlab = matlab.engine.start_matlab()
    curv_m = mlab.CenterlineCharacterization(path, filename, init_knots, order, nargout=1)
    #curv_m = mlab.CenterlineCharacterization(path, filename, init_knots, order, init_array_matlab, nargout=1)
    n = len(curv_m)
    curv_p = np.zeros(n)

    for i in range(n):
        curv_p[i] = curv_m[i][0]

    if plot: 
        plt.plot(curv_p, label="Knot free regression spline")
        plt.legend()
        plt.show()

    return curv_p



def discrete_geometry(line, neigh=10):
    len_line = line.GetNumberOfPoints()
    N = line.GetNumberOfPoints() 

    # Compute cumulative chord length
    t = np.zeros(N)
    p = []
    for i in range(N):
        p.append(np.array(list(line.GetPoint(i))))
        p[i] = np.array(p[i])

    norms = [la.norm(p[j] - p[j-1]) for j in range(1,N)]
    s = sum(norms)
    for i in range(1, N):
        s1 = sum(norms[:i+1])
        t[i] = s1 / s

    # Radius of sliding neighbourhood
    m = neigh

    dxdt = np.zeros(N)
    dydt = np.zeros(N)
    dzdt = np.zeros(N)

    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)

    for i in range(N):
        x[i] = p[i][0]
        y[i] = p[i][1]
        z[i] = p[i][2]
 

    for i in range(0, m):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(0, 2*m+1)])
        dxdt[i] = sum([ (t[j] - t[i])*(x[j]-x[i]) for j in range(0, 2*m+1)]) / t_sum
        dydt[i] = sum([ (t[j] - t[i])*(y[j]-y[i]) for j in range(0, 2*m+1)]) / t_sum
        dzdt[i] = sum([ (t[j] - t[i])*(z[j]-z[i]) for j in range(0, 2*m+1)]) / t_sum

    for i in range(m, N-m):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(i-m, i+m+1)])
        dxdt[i] = sum([ (t[j] - t[i])*(x[j]-x[i]) for j in range(i-m, i+m+1)]) / t_sum
        dydt[i] = sum([ (t[j] - t[i])*(y[j]-y[i]) for j in range(i-m, i+m+1)]) / t_sum
        dzdt[i] = sum([ (t[j] - t[i])*(z[j]-z[i]) for j in range(i-m, i+m+1)]) / t_sum

    for i in range(N-m, N):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(N-2*m, N)])
        dxdt[i] = sum([ (t[j] - t[i])*(x[j]-x[i]) for j in range(N-2*m-1, N)]) / t_sum
        dydt[i] = sum([ (t[j] - t[i])*(y[j]-y[i]) for j in range(N-2*m-1, N)]) / t_sum
        dzdt[i] = sum([ (t[j] - t[i])*(z[j]-z[i]) for j in range(N-2*m-1, N)]) / t_sum

    dgammadt = []
    dgammadt_norm = np.zeros(N)
    for i in range(N):
        dgammadt.append(np.array([dxdt[i], dydt[i], dzdt[i]]))
        dgammadt_norm[i] = la.norm(dgammadt[i])
    
    tg = []
    for i in range(N):
        tg.append(dgammadt[i] / dgammadt_norm[i])

    t1 = np.zeros(N)
    t2 = np.zeros(N)
    t3 = np.zeros(N)

    for i in range(N):
        t1[i] = tg[i][0]
        t2[i] = tg[i][1]
        t3[i] = tg[i][2]
    
    dt1dt = np.zeros(N)
    dt2dt = np.zeros(N)
    dt3dt = np.zeros(N)
    #m #= m+1


    for i in range(0, m):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(0, 2*m+1)])
        dt1dt[i] = sum([ (t[j] - t[i])*(t1[j]-t1[i]) for j in range(0, 2*m+1)]) / t_sum
        dt2dt[i] = sum([ (t[j] - t[i])*(t2[j]-t2[i]) for j in range(0, 2*m+1)]) / t_sum
        dt3dt[i] = sum([ (t[j] - t[i])*(t3[j]-t3[i]) for j in range(0, 2*m+1)]) / t_sum

    for i in range(m, N-m):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(i-m, i+m+1)])
        dt1dt[i] = sum([ (t[j] - t[i])*(t1[j]-t1[i]) for j in range(i-m, i+m+1)]) / t_sum
        dt2dt[i] = sum([ (t[j] - t[i])*(t2[j]-t2[i]) for j in range(i-m, i+m+1)]) / t_sum
        dt3dt[i] = sum([ (t[j] - t[i])*(t3[j]-t3[i]) for j in range(i-m, i+m+1)]) / t_sum
    
    for i in range(N-m, N):
        t_sum = sum([ (t[j] - t[i])**2 for j in range(N-2*m, N)])
        dt1dt[i] = sum([ (t[j] - t[i])*(t1[j]-t1[i]) for j in range(N-2*m-1, N)]) / t_sum
        dt2dt[i] = sum([ (t[j] - t[i])*(t2[j]-t2[i]) for j in range(N-2*m-1, N)]) / t_sum
        dt3dt[i] = sum([ (t[j] - t[i])*(t3[j]-t3[i]) for j in range(N-2*m-1, N)]) / t_sum

    dtgdt = []
    dtgdt_norm = np.zeros(N)
    for i in range(N):
        dtgdt.append(np.array([dt1dt[i], dt2dt[i], dt3dt[i]]))
        dtgdt_norm[i] = la.norm(dtgdt[i])

    curv = np.zeros(N)
    for i in range(N):
        curv[i] = dtgdt_norm[i] / dgammadt_norm[i]

    curv = resample(curv, len_line)

    return line, curv


def writeCenterline(centerline):
    P = []
    for i in range(centerline.GetNumberOfPoints()):
        P.append(np.array(centerline.GetPoint(i)))

    clfile = "centerline.txt"
    with open(clfile,'wb') as f:
        for p in P:
            f.write("%s %s %s\n" % (p[0],p[1],p[2]))
    return clfile

def landmarking_bogunovic(centerline, folder, showplot, resamp_step, smooth=False, nknots=25):
    """Takes the file path to a ceterline patch created by
    patchandinterpolatecenterlines.CreateParentArteryPatches, and spline the
    centerline and uses Bogunevic et al. (2012) to do automated landmarking"""
    print "Case: %s" % folder
    vmtk = False
    knotfree = False
    spline = True
    resamp = False
    disc = False
    if knotfree:
        init = nknots
        order = 5

        # Get attributes 
        line = CenterlineAttribiutes(centerline, smooth=False)
        clfile = writeCenterline(line)
        curvature_ = spline_matlab(".", clfile, init_knots=init, order=float(5),plot=False)
        curvature__ = gauss(curvature_,2)
        curvature = []
        for c in curvature__:
            curvature.append([c])
        curvature = np.array(curvature)

    elif vmtk: 
        neigh = 20
        #line = CenterlineResampling(centerline, 0.1)
        #line = CenterlineGeometry(line, outputsmoothed=1, factor=1.8, iterations=200)
        line = centerline
        line = CenterlineAttribiutes(line, smooth=True, factor=0.5, it=100)
        line, curvature__ = discrete_geometry( line, neigh=neigh)
        curvature = []
        for c in curvature__:
            curvature.append([c])
        curvature = np.array(curvature)

    elif disc:
        neigh = 20
        line = CenterlineAttribiutes(centerline, smooth=False)
        line, curvature__ = discrete_geometry(line, neigh=neigh)
        curvature = []
        for c in curvature__:
            curvature.append([c])
        curvature = np.array(curvature)

    elif spline: 
        nknots = 11
        if resamp:
            resamp_step = 0.1
            centerline = CenterlineResampling(centerline, length=resamp_step)

        line, max_point_ids, min_point_ids = splineCenterline(centerline, smooth, nknots)
       # line = CenterlineAttribiutes(line, smooth=False)
        torsion = get_array("Torsion", line)
        torsion = gaussian_filter(torsion , 5)

        curvature = get_array("Curvature", line)

    if not spline:
        # Make k1-k2 basis
        max_point_ids = list(argrelextrema(curvature, np.greater)[0])
        min_point_ids = list(argrelextrema(curvature, np.less)[0])

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
        N = get_array("FrenetNormal", line, k=3)
        #curvature = get_array("Curvature", line)

        k2 = (curvature.T * (E1[:,1]*N[:,0] - N[:,1]*E1[:,0]) / \
                            (E2[:,1]*E1[:,0] - E2[:,0]*E1[:,1]))[0]
        k1 = (-(curvature.T * N[:,0] + k2*E2[:,0]) / E1[:,0])[0]

        for k in [(k1, "k1"), (k2, "k2")]:
            k_array = create_vtk_array(k[0], k[1])
            line.GetPointData().AddArray(k_array)


    length = get_curvilinear_coordinate(line)
    k1 = get_array("k1", line)
    k2 = get_array("k2", line)

    # Remove a min / max point that is in reality a saddle point
    for i in min_point_ids:
        for j in max_point_ids:
            if abs(i-j) < 5 and abs(curvature[i] - curvature[j]) < 0.01:
                if i in min_point_ids:
                    min_point_ids.remove(i)
                if j in max_point_ids:
                    max_point_ids.remove(j)

    #from IPython import embed; embed()
    k1_points = k1[max_point_ids]
    k2_points = k2[max_point_ids]
    k_points = np.zeros((k1_points.shape[0], 2))
    k_points[:,0] = k1_points[:,0]
    k_points[:,1] = k2_points[:,0]
    tetha = np.zeros(k1_points.shape[0] - 1)
    for i in range(tetha.shape[0]):
        a = k_points[i,:] / np.sqrt(np.sum(k_points[i,:]**2))
        b = k_points[i+1,:] / np.sqrt(np.sum(k_points[i+1,:]**2))
        tetha[i] = math.acos(np.dot(a, b))
        tetha[i] = tetha[i] * 180 / math.pi
    
    Z = np.zeros(length.shape[0])
    Y = np.zeros(length.shape[0])
    X = np.zeros(length.shape[0])
    for i in range(Z.shape[0]):
        X[i] = line.GetPoints().GetPoint(i)[0]
        Y[i] = line.GetPoints().GetPoint(i)[1]
        Z[i] = line.GetPoints().GetPoint(i)[2]

    # Tolerance parameters from Bogunevic et al. (2012)
    tol_ant_post = 60
    tol_sup_ant = 45
    tol_post_inf = 45
    tol_inf_end = 110

    #from matplotlib.pyplot import plot, hold, show, xlabel, ylabel, legend, axis
    #plot(length, Z) #[length[m], length[m]], [0,1])
    #show()
    
    # Find max coronal coordinate
    value_index = Z[argrelextrema(Z, np.less)[0]].min()
    max_coronal_ids = np.array(Z.tolist().index(value_index))
    if abs(length[max_coronal_ids]-length[-1]) > 30:
        print "Sanity check failed"
        #sys.exit(0)
        return None

    # Plot centerline

    """
    hold("on")
    plot(length[max_point_ids], curvature[max_point_ids], 'g^', label="Max curv")
    plot(length[min_point_ids], curvature[min_point_ids], 'b^', label="Min curv")
    legend()
    show()
    """

    def find_interface(start, dir, tol, part):
        stop = dir if dir == -1 else tetha.shape[0]
        sucess = False
        for i in range(start-1, stop, dir):
            if tetha[i] > tol:
                sucess = True
                break

        if sucess:
            start = max_point_ids[i]
            stop = max_point_ids[i + 1]
            b = min_point_ids[i+1]
            index = ((min_point_ids > start) * (min_point_ids < stop)).nonzero()[0]

            min_point = min_point_ids[index]
            interfaces[part] = min_point

        elif not sucess and part == "sup_ant":
            print "Where not able to identify the interface between the" + \
                  "anterior and superior bend. Chekc the coronal coordinates"
            return None

        elif not sucess and part != "inf_end":
            print "The geometry is to short to be classified with superior" + \
                    ", anterior, posterior and inferior."
            return None

        elif not sucess and part == "inf_end":
            interfaces["inf_end"] = np.array([0])
            i = 0
            print "End of inferior is at the end of the geometry, this might" + \
                  "affect the geometry stats"
        else:
            print "Something happend, idea: some bend ended at the last point"
            print part, min_point
            return None

        return i

    interfaces = {}
    min_point_ids = np.array(min_point_ids)
    #min_point_ids = min_point_ids[::-1]
    #max_point_ids = max_point_ids[::-1]
    #print min_point_ids
    index = np.array((max_coronal_ids > max_point_ids).nonzero()[0]).max()
    start = find_interface(index, -1, tol_ant_post, "ant_post")
    if start is None:
        return None
    start = find_interface(start, -1, tol_post_inf, "post_inf")
    if start is None:
        return None
    start = find_interface(start, -1, tol_inf_end, "inf_end")
    start = find_interface(index + 1, 1, tol_sup_ant, "sup_ant")
    if start is None:
        return None

    # Find a "center" of each bend
    bends = ["inferior", "posterior", "anterior", "superior"]
    values = [interfaces["inf_end"], interfaces["post_inf"],
              interfaces["ant_post"], interfaces["sup_ant"], 
              np.array([curvature.shape[0]])]

    max_landmarks = {}
    for i in range(len(bends)):
        """
        curv_part = curvature[values[i][0]: values[i+1][0] + 1][:, 0]
        max_new = []

        for j in range(len(max_point_ids)):
            if values[i] < max_point_ids[j] < values[i+1] + 1:
                max_new.append(max_point_ids[j])
        max_new = np.array(max_new)

        while max_new.shape[0] > 1:
            Gauss = gaussian(curv_part.shape[0], std=curv_part.shape[0]//2)    
            new = np.convolve(curv_part, Gauss, 'same')
            max_new = argrelextrema(new, np.greater)[0]

        #max_landmarks[bends[i]] = max_new + values[i]
        """

        # Following Bogunovic
        #max_landmarks.pop("superior")
        landmarks = {}
        for k, v in max_landmarks.iteritems():
            landmarks[k] = line.GetPoints().GetPoint(int(v))
        for k, v in interfaces.iteritems():
            landmarks[k] = line.GetPoints().GetPoint(int(v))

    print "Case = %s, Length = %.5f" % (folder, length[-1])
    #return landmarks
    method = "bogunovic"
    smoothing = "smooth" if smooth else "spline"
    nknots = resamp_step if smooth else nknots
    if landmarks is not None:
        WritePolyData(line, path.join(folder, "%s_%s_%s_resamp01_11.vtp" % (method, nknots, smoothing)))
        writeParameters(landmarks, folder)
        createParticles(folder, method, nknots, smoothing, resamp_step)



def createParticles(folder, method, nknots, smooth, resamp_step):
    # Create a file with points where bends are located and 
    # remove points from manifest

    manifest = folder + "/manifest.txt"
    filename = folder + "/bogbog.particles" 

    output = open(filename, "w")
    mani= open(manifest, "r")
    lines = mani.readlines()
    mani.close()
    mani = open(manifest, "w")

    for line in lines:
        if line.split()[1][0] == "(":
            coord = line.split()[1:]
            point = "%s %s %s" % (coord[0][1:-1], coord[1][:-1], coord[2][:-1])
            #print "Location: %s - " % line.split()[:1],point
            output.write(point + "\n")
        else:
            mani.write(line)

    output.close()
    mani.close()

endids = []

def main(folder):
    methods=["bogunovic"]#, "piccinelli"]
    #methods=[ "piccinelli"]
    smoothing=[False]
    showplot=False
    nknots_val_bog = [11]
    case = folder.split("/")[-1]
    nknots_val_pic =[]# [5,11,25]
    resamp_steps = [ 0.2, 0.3,0.5]#,0.3, 0.5]

    centerline_path = path.join(folder, "surface", "model_usr_centerline.vtp")
    centerline_bif_path = path.join(folder, "surface", "model_usr_centerline_bif.vtp")
    centerline_bif = ReadPolyData(centerline_bif_path)
    centerline = ReadPolyData(centerline_path)

    centerlineSpacing = math.sqrt(vtk.vtkMath.Distance2BetweenPoints( \
                                  centerline.GetPoint(2), \
                                  centerline.GetPoint(3)))
    divergingTolerance = centerlineSpacing / divergingRatioToSpacingTolerance

    data = getData(centerline, centerline_bif, centerlineSpacing)
    endid = data["bif"]["ID_div"]
    endids.append(endid)
    line = ExtractSingleLine(centerline, 0, startID=0, endID=data["bif"]["ID_div"])

    WritePolyData(line, path.join(folder, "surface", "carotid_siphon.vtp"))
    fold = folder[30:]
    for smooth in smoothing:
        for method in methods:
            if method == "bogunovic": 
                if smooth:
                    for resamp_step in resamp_steps:
                        if fold in ["C0037","P0220", "P0251"]:
                            nknots = 12
                        elif fold in ["P0250"]:
                            nknots = 14
                        elif fold in ["C0015"]:
                            nknots = 10
                        else:
                            nknots = 11
              #          if fold in ["P0220"]:
                        smooth=False
                        landmarking_bogunovic(line, folder,showplot, resamp_step, smooth ,nknots=nknots)
                else: 
                    for nknots in nknots_val_bog:
                        resamp_step = 0.3
                        if fold in ["P0251"]:
                            nknots = 12
                        elif fold in ["P0250"]:
                            nknots = 14
                        elif fold in ["P0086"]:
                            nknots = 10
                        
                        landmarking_bogunovic(line, folder,showplot, resamp_step, smooth ,nknots)

            if method == "piccinelli":
                if smooth:
     #               for resamp_step in resamp_steps:
                    if fold in ["P0157", "P0207", "P0220"]:
                        resamp_step = 0.3
                    else:
                        resamp_step = 0.2
                    landmarking_piccinelli(line, folder,showplot, resamp_step, smooth ,nknots=20)
                else: 
         #           for nknots in nknots_val_pic:
                    resamp_step = 0.2
                    landmarking_piccinelli(line, folder,showplot, resamp_step, smooth ,nknots)

    
    #landmarks = landmarking(line, folder, method)
    #if landmarks is not None:
    #    writeParameters(landmarks, folder)
    

def check_landmarked(dirpath):
    parameters = getParameters(dirpath)
    return parameters.has_key("sup_ant")


if __name__ == '__main__':
    maindir = '/home/henrik/master_uio/cases/'
    for folder in sorted(listdir(maindir)):
        dirpath = path.join(maindir, folder)
        if path.isdir(dirpath) and folder != "backup" and folder != ".git" and "test" not in folder:
            #if check_landmarked(dirpath):
            #    print "is allready marked, moving on!"
            #elif folder in ["C0075", "C0063", "P0134", "C0035", "C0086",
            #        "C0088a", "C0067", "C0085", "P0252", "C0005", "C0022", "C0074a"]:
            #    print "special case, moving on!"
            #else:
            #    print "starting to landmark..."
            #if folder in ["p0252"]:#["C0015", "p0228", "p0250"]:
            main(dirpath)

    #with open('end_ids.txt', 'a') as the_file:
    #    for ID in endids:
    #        the_file.write('%i\n' % ID)

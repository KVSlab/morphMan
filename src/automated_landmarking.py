# -*- coding: utf-8 -*-
#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import sys
import math

from argparse import ArgumentParser
from common import *
from patchandinterpolatecenterlines import *
from os import path, listdir
from scipy.interpolate import splrep, splev
from scipy.ndimage.filters import gaussian_filter
from subprocess import check_output, STDOUT
from IPython import embed

def read_command_line():
    """
    Read arguments from commandline
    """
    parser = ArgumentParser()

    parser.add_argument('-d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('-c','--case', type=str, default=None, help="Choose case")
    parser.add_argument('-cm','--curv_method', type=str, default="spline",
            help="Choose which method used for computing curvature: spline (default) | vmtk | disc |")

    parser.add_argument('-a','--algorithm', type=str, default="bogunovic",
            help="Choose which landmarking algorithm to use: bogunovic or piccinelli")
    parser.add_argument('-r','--resamp_step', type=float, default=None,
                        help="Choose if the centerline should be resampled with this resampling step")
    parser.add_argument('-nk','--nknots', type=int, default=11,
                        help="Number of knots used in B-splines.")
    parser.add_argument('-s', '--smooth', type=bool, default=False,
                        help="If the original centerline should be smoothed " +\
                        "when computing the centerline attribiutes", metavar="smooth")
    parser.add_argument('-facc','--factor_curvature', type=float, default=1.0,
                        help="Smoothing factor for computing curvature.")
    parser.add_argument('-fact','--factor_torsion', type=float, default=1.0,
                        help="Smoothing factor for computing torsion.")
    parser.add_argument('-it','--iterations', type=int, default=100,
                        help="Smoothing iterations.")
    args = parser.parse_args()

    return args.dir_path, args.case, args.curv_method, args.resamp_step, args.smooth, args.algorithm, args.nknots, args.factor_curvature, args.factor_torsion, args.iterations

def landmarking_bogunovic(centerline, folder, curv_method, algorithm,
        resamp_step, smooth=False, nknots=25, smoothingfactor=1.0, iterations=100):
    """
    Perform landmarking of an input centerline to
    identify different segments along the vessel.
    Landmarking algorithm based on Bogunevic et al. (2012).
    Bends are identified where the angle in the k1-k2 basis
    exceeds a certain limit.

    Args:
        centerline (vtkPolyData): Centerline data points.
        folder (str): Location of case to landmark.
        curv_method (str): Method used for computing curvature.
        algorithm (str): Name of landmarking algorithm.
        resamp_step (float): Resampling step. Is None if no resampling.
        smooth (bool): Smoothes centerline with VMTK if True.
        nknots (int): Number of knots for B-splines.
        smoothingfactor (float): Smoothing factor used in VMTK
        iterations (int): Number of smoothing iterations.
    """

    if resamp_step is not None:
        centerline = vmtk_centerline_resampling(centerline, length=resamp_step)

    elif curv_method == "vmtk":
        line = vmtk_centerline_attributes(centerline)
        line = vmtk_centerline_geometry(line,smooth, factor=smoothingfactor, iterations=iterations)
        curvature_ = get_array("Curvature", line)
        curvature__ = gaussian_filter(curvature_, 5)
        curvature = []
        for c in curvature__:
            curvature.append(c)
        curvature = np.array(curvature)

    elif curv_method == "disc":
        neigh = 20
        line = vmtk_centerline_attributes(centerline)
        line = vmtk_centerline_geometry(line,smooth, factor=smoothingfactor, iterations=iterations)
        line, curvature__ = discrete_geometry(line, neigh=neigh)
        curvature = []
        for c in curvature__:
            curvature.append([c])
        curvature = np.array(curvature)

    elif curv_method == "spline":
        line, max_point_ids, min_point_ids = spline_and_geometry(centerline, smooth, nknots)
        curvature = get_array("Curvature", line)

    if curv_method != "spline":
        max_point_ids = list(argrelextrema(curvature, np.greater)[0])
        min_point_ids = list(argrelextrema(curvature, np.less)[0])
        get_k1k2_basis(curvature, line)

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

    # Find max coronal coordinate
    value_index = Z[argrelextrema(Z, np.less)[0]].min()
    max_coronal_ids = np.array(Z.tolist().index(value_index))
    if abs(length[max_coronal_ids]-length[-1]) > 30:
        print("Sanity check failed")
        return None

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
            #min_point = min(min_point_ids[index])
            min_point = (min_point_ids[index])
            interfaces[part] = min_point

        elif not sucess and part == "sup_ant":
            print("Where not able to identify the interface between the" + \
                  "anterior and superior bend. Chekc the coronal coordinates")
            return None

        elif not sucess and part != "inf_end":
            print("The geometry is to short to be classified with superior" + \
                    ", anterior, posterior and inferior.")
            return None

        elif not sucess and part == "inf_end":
            interfaces["inf_end"] = np.array([0])
            i = 0
            print("End of inferior is at the end of the geometry, this might" + \
                  "affect the geometry stats")
        else:
            print("Something happend, idea: some bend ended at the last point")
            print(part, min_point)
            return None

        return i

    interfaces = {}
    min_point_ids = np.array(min_point_ids)

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
        # Following Bogunovic
        landmarks = {}
        for k, v in max_landmarks.iteritems():
            landmarks[k] = line.GetPoints().GetPoint(int(v))
        for k, v in interfaces.iteritems():
            landmarks[k] = line.GetPoints().GetPoint(int(v))


    # Save landmarks
    print("Case %s was successfully landmarked." % (folder[-5:]))
    if landmarks is not None:
        write_parameters(landmarks, folder)
        create_particles(folder, algorithm, curv_method)


def landmarking_piccinelli(centerline, folder, curv_method, algorithm, resamp_step=None,
                           smooth=False, nknots=None, factor_curv=None, factor_tor=None,
                           iterations=100);
    """
    Perform landmarking of an input centerline to
    identify different segments along the vessel.
    Landmarking algorithm based on Bogunevic et al. (2012).
    Uses curvature and torsion to objectivley subdivide
    the siphon into bends.
    Subdivision of individual siphon bends is
    performed by identifying locations of curvature and torsion peaks along
    the siphon and defining a bend for each curvature peak delimited by 2
    enclosing (proximal and distal) torsion peaks.

    Args:
        centerline (vtkPolyData): Centerline data points.
        folder (str): Location of case to landmark.
        curv_method (str): Method used for computing curvature.
        algorithm (str): Name of landmarking algorithm.
        resamp_step (float): Resampling step. Is None if no resampling.
        smooth (bool): Smoothes centerline with VMTK if True.
        nknots (int): Number of knots for B-splines.
        factor_curv (float): Smoothing factor for computing curvature.
        factor_torsion (float): Smoothing factor for computing torsion.
        iterations (int): Number of smoothing iterations.
    """
    print("Case: %s" % folder.split("/")[-1])

    if resamp_step != None:
        centerline = CenterlineResampling(centerline, step)

    if curv_method == "spline":
        line, max_point_ids, min_point_ids = spline_and_geometry(centerline, smooth, nknots)

        # Get curvature and torsion, find peaks
        curvature = get_array("Curvature", line)
        torsion = get_array("Torsion", line)
        torsion_smooth = gaussian_filter(torsion , 10)
        max_point_tor_ids = list(argrelextrema(abs(torsion_smooth), np.greater)[0])

    elif curv_method == "vmtk":
        line = vmtk_centerline_geometry(centerline, True, outputsmoothed=0,
                                        factor=factor_curv, iterations=iterations)
        line_tor = vmtk_centerline_geometry(centerline, True, outputsmoothed=0,
                                            factor=factor_tor, iterations=iterations)

        # Get curvature and torsion, find peaks
        curvature = get_array("Curvature", line)
        torsion = get_array("Torsion", line_tor)
        torsion_smooth = gaussian_filter(torsion , 10)
        curvature_smooth = gaussian_filter(curvature, 10)
        max_point_ids = list(argrelextrema(curvature_smooth, np.greater)[0])
        max_point_tor_ids = list(argrelextrema(abs(torsion_smooth), np.greater)[0])

    else:
        print("Not a valid method for computing curvature/torsion")
        print("Methods available: 'spline' or 'vmtk'")
        sys.exit(1)

    # Extract local curvature minimums
    length = get_curvilinear_coordinate(line)

    # Remove points too close to the ends of the siphon
    for i in max_point_ids:
        if length[i] in length[-10:] or length[i] in length[:10]:
            max_point_ids.remove(i)

    # Remove curvature and torsion peaks too close to each other
    tolerance = 50
    dist = []
    dist_tor = []
    for i in range(len(max_point_ids)-1):
        dist.append(max_point_ids[i+1] - max_point_ids[i])
    for i in range(len(max_point_tor_ids)-1):
        dist_tor.append(max_point_tor_ids[i+1] - max_point_tor_ids[i])

    for i, dx in enumerate(dist):
        if dx < tolerance:
            curv1 = curvature[max_point_ids[i]]
            curv2 = curvature[max_point_ids[i+1]]
            if curv1 > curv2:
                max_point_ids[i+1] = None
            else:
                max_point_ids[i] = None

    for i, dx in enumerate(dist_tor):
        if dx < tolerance:
            tor1 = torsion_smooth[max_point_tor_ids[i]]
            tor2 = torsion_smooth[max_point_tor_ids[i+1]]
            if tor1 > tor2:
                max_point_tor_ids[i+1] = None
            else:
                max_point_tor_ids[i] = None

    max_point_ids = [ID for ID in max_point_ids if ID != None]
    max_point_tor_ids = [ID for ID in max_point_tor_ids if ID != None]

    # Define bend interfaces based on Piccinelli et al.
    def find_interface():
        interface = {}
        k=0
        start_id = 0
        for c in max_point_ids:
            for i in range(start_id, len(max_point_tor_ids) - 1):
                if max_point_tor_ids[i] < c and c < max_point_tor_ids[i+1] and not found:
                    interface["bend%s" % (k+1)] = np.array([max_point_tor_ids[i]])
                    k+=1
                    interface["bend%s" % (k+1)] = np.array([max_point_tor_ids[i+1]])
                    k+=1
                    start_id = i+1
                    found = True
            found = False

        return interface

    # Compute and extract interface points
    interfaces = find_interface()
    landmarks = {}
    for k, v in interfaces.iteritems():
        landmarks[k] = line.GetPoints().GetPoint(int(v))

    # Save landmarks
    if landmarks is not None:
        write_parameters(landmarks, folder)
        create_particles(folder, algorithm, curv_method)



def spline_and_geometry(line, smooth, nknots):
    """
    Compute attributes and geometric parameters of input
    centerline, using B-splines (SciPy).

    Args:
        line (vtkPolyData): Centerline data.
        smooth (bool): Smooth centerline with VMTK if True.
        nknots (int): Number of knots for B-splines.

    Returns:
        line (vtkPolyData): Splined centerline.
    Returns:
        max_point_ids (ndarray): Array of max curvature values
    Returns:
        min_point_ids (ndarray): Array of min curvature values
    """


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
    line = vmtk_centerline_attributes(line)
    line = vmtk_centerline_geometry(line, smooth)
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

    curvature = get_array("Curvature", line)
    line = get_k1k2_basis(curvature, line)

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



def create_particles(folder, algorithm, method):
    """
    Create a file with points where bends are located and
    remove points from manifest

    Args:
        folder (str): Case location.
        algorithm (str): Name of landmarking algorithm.
        method (str): Method used for computing curvature.
    """

    manifest = folder + "/info.txt"
    filename_all = folder + "/landmark_%s_%s.particles" % (algorithm, method)
    filename_siphon = folder + "/anterior_bend_%s_%s.particles"  % (algorithm, method)

    output_all = open(filename_all, "w")
    output_siphon = open(filename_siphon, "w")
    mani= open(manifest, "r")
    lines = mani.readlines()
    mani.close()

    for line in lines:
        if algorithm == "bogunovic":
            if line.split()[1][0] == "(":
                coord = line.split()[1:]
                point = "%s %s %s" % (coord[0][1:-1], coord[1][:-1], coord[2][:-1])
                output_all.write(point + "\n")
                if line.split()[0] == "sup_ant:" or line.split()[0] == "ant_post:":
                    output_siphon.write(point + "\n")

        elif algorithm == "piccinelli":
            if line.split()[0][:4] == "bend":
                coord = line.split()[1:]
                point = "%s %s %s" % (coord[0][1:-1], coord[1][:-1], coord[2][:-1])
                output_all.write(point + "\n")

    output_all.close()
    output_siphon.close()

def main(folder, curv_method, algorithm, resamp_step, nknots, smooth, factor_curv,
         factor_torsion, iterations):
    """
    Compute carotid siphon and perform landmarking.

    Args:
        folder (str): Location of case to landmark.
        curv_method (str): Method used for computing curvature.
        algorithm (str): Name of landmarking algorithm.
        resamp_step (float): Resampling step. Is None if no resampling.
        nknots (int): Number of knots for B-splines.
        smooth (bool): Smoothes centerline with VMTK if True.
        factor_curv (float): Smoothing factor used in VMTK for curvature
        factor_torsion (float): Smoothing factor used in VMTK for torsion
        iterations (int): Number of smoothing iterations.
    """
    case = folder.split("/")[-1]

    # Extract carotid siphon
    line = extract_carotid_siphon(folder)

    # Landmark
    if algorithm == "bogunovic":
        landmarking_bogunovic(line, folder, curv_method, algorithm, resamp_step, smooth ,nknots, factor_curv, iterations)
    elif algorithm == "piccinelli":
        landmarking_piccinelli(line, folder, curv_method, algorithm, resamp_step, smooth ,nknots, factor_curv, factor_torsion, iterations)
    else:
        print("Algorithm is not valid. Select between 'bogunovic' and 'piccinelli'.")



if __name__ == '__main__':
    basedir, case, curv_method, resamp_step, smooth, algorithm, nknots, factor_curv, factor_torsion, iterations = read_command_line()
    folders = listdir(basedir) if case is None else [case]
    for folder in sorted(folders):
        dirpath = path.join(basedir, folder)
        if path.isdir(dirpath):
            main(dirpath, curv_method, algorithm, resamp_step, nknots, smooth,
                    factor_curv, factor_torsion, iterations)

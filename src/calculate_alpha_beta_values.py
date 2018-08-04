# -*- coding: utf-8 -*-
from argparse import ArgumentParser
from os import listdir, path
from scipy.interpolate import griddata
from scipy import interpolate

import numpy as np
import numpy.linalg as la

def read_command_line():
    """
    Read arguments from commandline
    """
    parser = ArgumentParser()
    parser.add_argument('-d', '--dir_path', type=str, default=".",
                        help="Path to the folder with alpha-beta_vales")
    parser.add_argument("-p", "--param", type=str, default="curvature",
                        help="Parameter to compute.")

    args = parser.parse_args()
    return args.dir_path, args.param

def get_alpha_beta(dirpath, i, files, param):
    """
    Imports a matrix of parameter values corresponding
    to a (alpha,beta) point and perform spline interpolation
    to make a parameter surface.
    Parameter surface is used to find intersection with
    initial parameter value plus / minus one standard deviation.
    Three criterias to find suggested alpha-beta value from intersection.

    Args:
        dirpath (str): Location where parameter value data is stored.
        i (int): Index of case iteration
        files (str): Filename of parameter value data file.
        param (str): Parameter name.
    """

    files = path.join(dirpath, files)
    case = files[-9:-4]
    print "Working on case %s" % case

    # Get boundaries
    alphabeta_bound = np.loadtxt("alphabeta_bound.txt")
    amin, amax, bmin, bmax = alphabeta_bound[i][0], alphabeta_bound[i][1], alphabeta_bound[i][2], alphabeta_bound[i][3]

    # Set standard deviations used to find intersetcion
    if param == "curvature":
        sd_curv = 0.045
    elif param == "angle":
        sd_angle = 19.1

    # Defined SD planes for curvature
    # Tolerance added for adjusting SD
    # if there are no intersections found
    def cpsd(x, y, tol=0.0):
        return curv0 + sd_curv - tol

    def cmsd(x, y, tol=0.0):
        return curv0 - sd_curv + tol

    def curv_init(x, y, tol=0.0):
        return curv0

    def apsd(x, y, tol=0.0):
        return angle0 + sd_angle - tol

    def amsd(x, y, tol=0.0):
        return angle0 - sd_angle + tol

    def angle_init(x, y, tol=0.0):
        return angle0

    # Extract values
    with open(files) as file:
        data = [[float(digit) for digit in line.split()] for line in file]

    N = len(data)
    edge = [amin, amax, bmin, bmax]
    alpha = np.linspace(amin, amax, N)
    beta = np.linspace(bmin, bmax, N)
    alpha_long = np.linspace(amin, amax, 300)
    beta_long = np.linspace(bmin, bmax, 300)
    yy, xx = np.meshgrid(beta, alpha)

    points = np.zeros((N, 2))
    for i in range(len(xx)):
        points[i] = [alpha[i], beta[i]]

    # Spline interpolation
    f = interpolate.interp2d(beta, alpha, data, kind='cubic')
    if param == "curvature":
        curv0 = f(0, 0)
    elif param == "angle":
        angle0 = f(0, 0)
    zz = f(beta, alpha)

    Z = []
    if param == "curvature":
        methods = [cpsd, cmsd, curv_init]
    elif param == "angle":
        methods = [apsd, amsd, angle_init]

    # Find intersecting points
    # Reduces SD if no points are found
    for plane in methods:
        print "Method: %s" % (plane.__name__)

        zeros = alpha_beta_intersection(plane, f, alpha_long, beta_long)
        if len(zeros) > 10:
            dx = int(len(zeros) / 10.)
        elif len(zeros) > 0:
            dx = 1
        else:
            empty = True
            tol = 0.005 if param == "curvature" else 0.1
            maxiter = 50
            iterations = 0

            print "Found no points..Adjusting SD"
            while empty == True and iterations < maxiter:
                print "Iterations: %i" % (iterations + 1)
                zeros = alpha_beta_intersection(plane, f, alpha_long, beta_long, tol)
                if len(zeros) > 0:
                    empty = False
                iterations += 1
                if param == "curvature":
                    tol += 0.001
                elif param == "angle":
                    tol += 0.2

            if len(zeros) > 10:
                dx = int(len(zeros) / 10.)
            elif len(zeros) > 0:
                dx = 1

        # Check points and apply criterias
        # to find suggested values for alpha and beta
        if len(zeros) > 0:
            points = []
            for p in zeros:
                if plane.__name__ in ["cpsd", "amsd"]:
                    if p[1] < 0:
                        points.append(p)
                elif plane.__name__ in ["cmsd", "apsd"]:
                    if p[1] > 0:
                        points.append(p)

            if plane.__name__ not in ["curv_init", "angle_init"]:
                P = points[0]
                min_rad = 0.15
                dist = 10
                for p in points[1:]:
                    dist_tmp = la.norm(np.array(p))
                    if dist_tmp < dist and dist_tmp > min_rad:
                        dist = dist_tmp
                        P = p
                Z.append(P)

                # Write points to file
                write_alpha_beta_point(case, P, plane.__name__)


def write_alpha_beta_point(case, P, method):
    """
    Write suggested choice of (alpha, beta) to file.

    Args:
        case (str): Name of case.
        P (ndarray): Array containing alpha and beta value
        method (str): Info about parameter, and increase / decrease of parameter.
    """
    dirpath = path.join("alphabeta_values.txt")
    alpha = P[0]
    beta = P[1]
    with open(dirpath, "a") as f:
        f.write("%s %s alpha=%s beta=%s\n" % (case, method, alpha, beta))

def alpha_beta_intersection(method, f, alphas, betas, tol=0.0):
    """
    Iterate through values of alpha and beta
    and find intersection points with a given tolerance
    between initial value and initial value plus/minus SD.

    Args:
        method (function): Plane defining initial value plus/minus SD.
        f (interp2d): Interpolated surface of curvature or angle values.
        alphas (ndarray): List of alpha values.
        betas (ndarray): List of beta values.
        tol (float): Tolerance used for adjusting SD if needed.

    Returns:
        zeros (ndarray): Array containing (alpha,beta) tuples which intersect.
    """
    zeros = []
    for i in alphas:
        for j in betas:
            diff = abs(f(j, i) - method(0, 0, tol))
            if "c" in method.__name__:
                if diff < 0.001:
                    zeros.append([i, j])
            else:
                if diff < 0.05:
                    zeros.append([i, j])
    return zeros


def main(dirpath, param):
    """
    Get files containing parameter values
    and find suggested choice for alpha and beta,
    the compression / extension factors.

    Args:
        dirpath (str): Location of parameter text.
        param (str): Parameter to calculate.
    """
    files = listdir(dirpath)
    if param == "curvature":
        files_curv = sorted([f for f in listdir(dirpath) if f[4:8] in ["curv"]])
        for i in range(len(files_curv)):
            get_alpha_beta(dirpath, i, files_curv[i], param)
    elif param == "angle":
        files_angle = sorted([f for f in listdir(dirpath) if f[4:9] in ["angle"]])
        for i in range(len(files_angle)):
            get_alpha_beta(dirpath, i, files_angle[i], param)


if __name__ == "__main__":
    dirpath, param = read_command_line()
    main(dirpath, param)

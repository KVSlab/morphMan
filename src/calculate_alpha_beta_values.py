from argparse import RawDescriptionHelpFormatter, ArgumentParser
from scipy import interpolate

import numpy as np
import numpy.linalg as la
from common import get_path_names




def get_alpha_beta(input_filepath, quantity_to_compute, boundary, radius, limit):
    """
    Imports a matrix of parameter values corresponding
    to a (alpha,beta) point and perform spline interpolation
    to make a parameter surface.
    Parameter surface is used to find intersection with
    initial parameter value plus / minus one standard deviation.
    Three criterias to find suggested alpha-beta value from intersection.

    Args:
        input_filepath (str): Surface model filename and location where data is stored.
        quantity_to_compute(str): Parameter name.
        boundary (list): Boundary of searching grid.
        radius (float): Minimum radius of circle to search outside of.
        limit (float): Desired change in curvature / bend angle to achieve
    """

    # Get grid values
    base_path = get_path_names(input_filepath)
    grid_filepath = base_path + "_grid_values.txt"
    amin, amax, bmin, bmax = float(boundary[0]), float(boundary[1]), float(boundary[2]), float(boundary[3])

    # Set standard deviations used to find intersetcion
    if quantity_to_compute == "curvature":
        sd_curv = limit
    elif quantity_to_compute == "angle":
        sd_angle = limit

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
    with open(grid_filepath) as file:
        data = [[float(digit) for digit in line.split()] for line in file]

    N = len(data)
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
    if quantity_to_compute == "curvature":
        curv0 = f(0, 0)
    elif quantity_to_compute == "angle":
        angle0 = f(0, 0)
    methods = []
    if quantity_to_compute == "curvature":
        methods = [cpsd, cmsd, curv_init]
    elif quantity_to_compute == "angle":
        methods = [apsd, amsd, angle_init]

    # Find intersecting points
    # Reduces SD if no points are found
    for plane in methods:
        zeros = alpha_beta_intersection(plane, f, alpha_long, beta_long)
        if len(zeros) > 0:
            continue
        else:
            empty = True
            tol = 0.005 if quantity_to_compute == "curvature" else 0.1
            maxiter = 50
            iterations = 0

            print("Found no points..Adjusting SD")
            while empty and iterations < maxiter:
                print("Iterations: %i" % (iterations + 1))
                zeros = alpha_beta_intersection(plane, f, alpha_long, beta_long, tol)
                if len(zeros) > 0:
                    empty = False
                iterations += 1
                if quantity_to_compute == "curvature":
                    tol += 0.001
                elif quantity_to_compute == "angle":
                    tol += 0.2

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
                suggested_point = points[0]
                dist = 10
                for p in points[1:]:
                    dist_tmp = la.norm(np.array(p))
                    if radius < dist_tmp < dist:
                        dist = dist_tmp
                        suggested_point = p

                # Write points to file
                write_alpha_beta_point(base_path, suggested_point, plane.__name__)


def write_alpha_beta_point(base_path, suggested_point, method):
    """
    Write suggested choice of (alpha, beta) to file.

    Args:
        base_path (str): Path to file directory.
        suggested_point (ndarray): Array containing alpha and beta value
        method (str): Info about parameter, and increase / decrease of parameter.
    """
    save_path = base_path + "_alphabeta_values.txt"
    alpha = suggested_point[0]
    beta = suggested_point[1]
    with open(save_path, "a") as f:
        f.write("%s alpha=%s beta=%s\n" % (method, alpha, beta))


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



def read_command_line():
    """
    Read arguments from commandline
    """
    description = "Algorithm used to compute the recommended value for " + \
                  "alpha and beta based on surface interpolation and a " + \
                  "given limit, depening on the quantity to be computed. " + \
                  "Primarly implemented for computing angle and curvature values."

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)
    required = parser.add_argument_group('required named arguments')

    # Required arguments
    required.add_argument('-i', '--ifile', type=str, default=None,
                          help="Path to the surface model", required=True)
    required.add_argument("-q", "--quantity", type=str, default="curvature",
                          help="Parameter to compute. Choose between 'curvature' and 'angle'", required=True)
    # Optional arguments
    parser.add_argument('-r', '--radius', type=float, default=0.15,
                        help="Radius of bounding circle, limiting the choice of alpha and beta")
    parser.add_argument('-b', '--boundary', nargs='+', default=[-0.2, 1, -0.2, 1],
                        help='Bounds of grid, as a list: [alpha_min, alpha_max, beta_min, beta_max]')
    parser.add_argument('-l', '--limit', type=float, default=0.4,
                        help='Desired change in curvature / bend angle to achieve. Algorithm computes' +
                             'recommended values of alpha and beta for both plus and minus this change.')

    args = parser.parse_args()

    return dict(input_filepath=args.ifile, quantity_to_compute=args.quantity,
                radius=args.radius, boundary=args.boundary, limit=args.limit)


if __name__ == "__main__":
    get_alpha_beta(**read_command_line())

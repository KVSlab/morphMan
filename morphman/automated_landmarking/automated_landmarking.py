##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from morphman.automated_landmarking import *
from morphman.common import *


def automated_landmarking(input_filepath, approximation_method, resampling_step, algorithm, nknots, smooth_line,
                          smoothing_factor_curv, smoothing_factor_torsion, iterations, coronal_axis):
    """
    Compute ICA and perform landmarking of the input geometry.
    The following landmarking algorithms are available:
    - Piccinelli et al. (2011)
    - Bogunović et al. (2012)

    Args:
        input_filepath (str): Location of surface model to landmark.
        approximation_method (str): Method used for computing curvature and torsion.
        resampling_step (float): Resampling step. Is None if no resampling.
        algorithm (str): Name of landmarking algorithm.
        nknots (int): Number of knots for B-splines.
        smooth_line (bool): Smooths centerline with VMTK if True.
        smoothing_factor_curv (float): Smoothing factor used in VMTK for curvature
        smoothing_factor_torsion (float): Smoothing factor used in VMTK for torsion
        iterations (int): Number of smoothing iterations.
        coronal_axis (str) : Axis determining coronal coordinate (Bogunovic)
    """
    base_path = get_path_names(input_filepath)

    # Extract carotid siphon
    ica_centerline = extract_ica_centerline(base_path, input_filepath, resampling_step)

    # Check axial coordinate of centerline, reverse if needed
    ica_centerline = orient_centerline(ica_centerline)

    # Landmark
    if algorithm == "bogunovic":
        landmarking_bogunovic(ica_centerline, base_path, approximation_method, algorithm, resampling_step, smooth_line,
                              nknots, smoothing_factor_curv, iterations, coronal_axis)

    elif algorithm == "piccinelli":
        landmarking_piccinelli(ica_centerline, base_path, approximation_method, algorithm, resampling_step, smooth_line,
                               nknots, smoothing_factor_curv, smoothing_factor_torsion, iterations)


def read_command_line():
    """
    Read arguments from commandline
    """
    description = "Perform landmarking of an input centerline to" + \
                  "identify different segments along the vessel." + \
                  "Landmarking algorithm based on Bogunovic et al. (2012) " + \
                  " and Piccinelli et al. (2011)."

    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--ifile', type=str, default=None, required=True,
                          help="Path to the surface model")

    # Optional arguments
    parser.add_argument('-m', '--approximation-method', type=str, default="vmtk",
                        help="Choose which method used for computing curvature and torsion. Default is 'vmtk'.",
                        choices=['spline', 'vmtk', 'disc'])
    parser.add_argument('-a', '--algorithm', type=str, default="bogunovic",
                        help="Choose which landmarking algorithm to use: " +
                             "'bogunovic' or 'piccinelli'. Default is 'bogunovic'.",
                        choices=['bogunovic', 'piccinelli'])
    parser.add_argument('-ca', '--coronal-axis', type=str, default="z",
                        help="Axis describing coronal coordinate. Default is 'z'.",
                        choices=['x', 'y', 'z'])
    parser.add_argument('-k', '--nknots', type=int, default=11,
                        help="Number of knots used in B-splines.")
    parser.add_argument('-sl', '--smooth-line', type=str2bool, default=False,
                        help="If the original centerline should be smoothed " +
                             "when computing the centerline attributes")
    parser.add_argument('-facc', '--smoothing-factor-curvature', type=float, default=1.5,
                        help="Smoothing factor for computing curvature.")
    parser.add_argument('-fact', '--smoothing-factor-torsion', type=float, default=1.2,
                        help="Smoothing factor for computing torsion.")
    parser.add_argument('-it', '--iterations', type=int, default=100,
                        help="Smoothing iterations.")
    parser.add_argument("-r", "--resampling-step", type=float, default=0.1,
                        help="Resampling step in centerlines.")
    args = parser.parse_args()

    return dict(input_filepath=args.ifile, approximation_method=args.approximation_method,
                coronal_axis=args.coronal_axis, resampling_step=args.resampling_step, algorithm=args.algorithm,
                nknots=args.nknots, smooth_line=args.smooth_line, smoothing_factor_curv=args.smoothing_factor_curvature,
                smoothing_factor_torsion=args.smoothing_factor_torsion, iterations=args.iterations)


if __name__ == '__main__':
    automated_landmarking(**read_command_line())

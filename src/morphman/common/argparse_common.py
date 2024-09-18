##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.
## This software is distributed WITHOUT ANY WARRANTY; without even
## the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE. See the above copyright notices for more information.


from argparse import ArgumentTypeError


def str2bool(boolean):
    """Convert a string to boolean.

    Args:
        boolean (str): Input string.

    Returns:
        return (bool): Converted string.
    """
    if boolean.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif boolean.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise ValueError("Boolean value expected.")


def restricted_float(x):
    x = float(x)
    if x <= 0.0 or x >= 1.0:
        raise ArgumentTypeError("{} not in range [0.0, 1.0]".format(x))
    return x


def add_common_arguments(parser, required=True):
    # Required arguments
    if required:
        required = parser.add_argument_group("Required arguments")
    else:
        required = parser

    required.add_argument(
        "-i",
        "--ifile",
        type=str,
        default=None,
        required=True,
        help="Path to the surface model",
    )
    required.add_argument(
        "-o",
        "--ofile",
        type=str,
        default=None,
        required=True,
        help="Relative path to the output surface. The default folder is"
        + " the same as the input file, and a name with a combination of the"
        + " parameters.",
    )

    # General arguments
    parser.add_argument(
        "-s",
        "--smooth",
        type=str2bool,
        default=True,
        help="Smooth the voronoi diagram, default is False",
    )
    parser.add_argument(
        "-f",
        "--smooth-factor",
        type=float,
        default=0.25,
        help="If smooth option is true then each voronoi point"
        + " that has a radius less then MISR*(1-smooth_factor) at"
        + " the closest centerline point is removed.",
    )
    parser.add_argument(
        "-n",
        "--no-smooth",
        type=str2bool,
        default=False,
        help="If true and smooth is true the user, if no-smooth-point is"
        + " not given, the user can provide points where the surface not will"
        + " be smoothed.",
    )
    parser.add_argument(
        "--no-smooth-point",
        nargs="+",
        type=float,
        default=None,
        help="If model is smoothed the user can manually select points on"
        + " the surface that will not be smoothed. A centerline will be"
        + " created to the extra point, and the section were the centerline"
        + " differ from the other centerlines will be kept un-smoothed. This"
        + " can be particle for instance when manipulating geometries"
        + " with aneurysms",
    )
    parser.add_argument(
        "-b",
        "--poly-ball-size",
        nargs=3,
        type=int,
        default=[120, 120, 120],
        help="The size of the poly balls that will envelope the new"
        + " surface. The default value is 120, 120, 120. If two tubular"
        + " structures are very close compared to the bounds, the poly ball"
        + " size should be adjusted. For quick proto typing we"
        + " recommend ~100 in all directions, but >250 for a final "
        + " surface.",
        metavar="size",
    )
    parser.add_argument(
        "--resampling-step",
        type=float,
        default=0.1,
        help="Resampling step in centerlines.",
    )

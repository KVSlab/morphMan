##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from common import str2bool

def add_common_arguments(parser):
    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--ifile', type=str, default=None, required=True,
                          help="Path to the surface model")
    required.add_argument("-o", "--ofile", type=str, default=None, required=True,
                          help="Relative path to the output surface. The default folder is" +
                               " the same as the input file, and a name with a combination of the" +
                               " parameters.")

    # General arguments
    parser.add_argument("-m", "--method", type=str, default="variation",
                        choices=["variation", "stenosis", "area"],
                        help="Methods for manipulating the area in the region of interest:" +
                             "\n1) 'variation' will increase or decrease the changes in area" +
                             " along the centerline of the region of interest." +
                             "\n2) 'stenosis' will create or remove a local narrowing of the" +
                             " surface. If two points is provided, the area between these" +
                             " two points will be linearly interpolated to remove the narrowing." +
                             " If only one point is provided it is assumed to be the center of" +
                             " the stenosis. The new stenosis will have a sin shape, however, any" +
                             " other shape may be easly implemented." +
                             "\n3) 'area' will inflate or deflate the area in the region of" +
                             " interest.")
    parser.add_argument('-s', '--smooth', type=str2bool, default=True,
                        help="Smooth the voronoi diagram, default is False")
    parser.add_argument('-f', '--smooth_factor', type=float, default=0.25,
                        help="If smooth option is true then each voronoi point" +
                             " that has a radius less then MISR*(1-smooth_factor) at" +
                             " the closest centerline point is removed.")
    parser.add_argument("-n", "--no-smooth", type=bool, default=False,
                        help="If true and smooth is true the user, if no-smooth-point is" +
                             " not given, the user can provide points where the surface not will" +
                             " be smoothed.")
    parser.add_argument("--no-smooth-point", nargs="+", type=float, default=None,
                        help="If model is smoothed the user can manually select points on" +
                             " the surface that will not be smoothed. A centerline will be" +
                             " created to the extra point, and the section were the centerline" +
                             " differ from the other centerlines will be keept un-smoothed. This" +
                             " can be practicle for instance when manipulating geometries" +
                             " with aneurysms")
    parser.add_argument("-b", "--poly-ball-size", nargs=3, type=int, default=[120, 120, 120],
                        help="The size of the poly balls that will envelope the new" +
                             " surface. The default value is 120, 120, 120. If two tubular" +
                             " structures are very close compared to the bounds, the poly ball" +
                             " size should be adjusted. For quick proto typing we" +
                             " recommend ~100 in all directions, but >250 for a final " +
                             " surface.", metavar="size")
    parser.add_argument("--resampling-step", type=float, default=0.1,
                        help="Resampling step in centerlines.")

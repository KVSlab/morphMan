##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from common.argparse_common import *
from common.surface_operations import *
from common.vessel_reconstruction_tools import *


def manipulate_surface(input_filepath, output_filepath, smooth, smooth_factor, no_smooth,
                       no_smooth_point, poly_ball_size, resampling_step,
                       region_of_interest, region_points, upper, frequency,
                       frequency_deviation, noise, radius_max, absolute):
    """
    Controll the surface roughness by removing or adding points from the Voronoi diagram
    of the surface. A less general version of the smoothing algorithm was first presented
    in Ford et al. (2009).
     Args:
        input_filepath (str): Path to input surface.
        output_filepath (str): Path to output surface.
        smooth (bool): Determine if the voronoi diagram should be smoothed.
        smooth_factor (float): Smoothing factor used for voronoi diagram smoothing.
        resampling_step (float): Resampling step used to resample centerlines.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' ]
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        poly_ball_size (list): Resolution of polyballs used to create surface.
        upper (float): Upper smoothing factor
        frequency (float): Frequency at which noise is added to the voronoi diagram, based on points along the centerline
        frequency_deviation (float): Standard deviation of frequency
        absolute (bool): Absolute value for the smoothing criteria
        radius_max (float): Used to pick MISR multiplier to create noise on surface
        noise (bool): Turn on/off adding noise to the surface.
    """
    # Filenames
    base_path = get_path_names(input_filepath)

    # Output filepaths
    # Centerliens
    centerlines_path = base_path + "_centerline.vtp"

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Get inlet and outlets
    inlet, outlets = get_inlet_and_outlet_centers(surface, base_path)
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=resampling_step,
                                                         smooth=False, base_path=base_path)

    centerline_splined, centerline_remaining, centerline_diverging, region_points, diverging_ids = get_line_to_change(
        capped_surface, centerlines, region_of_interest, "bend", region_points, 0)

    # Split the Voronoi diagram
    if centerline_diverging is not None:
        centerline_remaining = vtk_merge_polydata([centerline_diverging] + centerline_remaining)

    centerline_regions = [centerline_splined, centerline_remaining]
    voronoi_relevant, voronoi_rest = get_split_voronoi_diagram(voronoi, centerline_regions)

    # Compute and smooth voronoi diagram (not aneurysm)
    if smooth:
        print("-- Smooth the Voronoi diagram.")
        if no_smooth:
            no_smooth_cl = get_no_smooth_cl(capped_surface, centerlines, base_path, smooth,
                                            no_smooth, voronoi, no_smooth_point, pole_ids,
                                            resampling_step, region_points)
        else:
            no_smooth_cl = None

        voronoi = smooth_voronoi_diagram(voronoi, centerlines, smooth_factor, no_smooth_cl,
                                         absolute, upper)

    if noise:
        print("-- Add noise to the Voronoi diagram.")
        voronoi = add_noise_to_voronoi_diagram(voronoi_relevant, centerline_regions, radius_max, frequency,
                                               frequency_deviation)

    # Write a new surface from the new voronoi diagram
    print("-- Create new surface.")
    new_surface = create_new_surface(voronoi, poly_ball_size)
    write_polydata(new_surface, "unprepared_surface.vtp")  # DEBUG

    print("-- Preparing surface for output.")
    new_surface = prepare_output_surface(new_surface, surface, centerlines,
                                         output_filepath, test_merge=False, changed=True,
                                         old_centerline=centerlines)

    print("-- Writing new surface to {}.".format(output_filepath))
    write_polydata(new_surface, output_filepath)


def add_noise_to_voronoi_diagram(voronoi, centerline, radius_max, frequency, frequency_deviation):
    """

    Args:
        voronoi (vtkPolyData): Voronoi Diagram to be smoothed
        centerline (vtkPolyData): Centerline along relevant Voronoi diagram
        frequency (float): Frequency at which noise is added to the voronoi diagram, based on points along the centerline
        frequency_deviation (float): Standard deviation of frequency
        radius_max (float): Used to pick MISR multiplier to create noise on surface

    Returns:
        vtkPolyData: Noisy Voronoi diagram
    """
    n = voronoi.GetNumberOfPoints()
    radius_array = get_vtk_array(radiusArrayName, 1, n)
    radius_array_data = voronoi.GetPointData().GetArray(radiusArrayName).GetTuple1

    new_voronoi = vtk.vtkPolyData()
    cell_array = vtk.vtkCellArray()
    points = vtk.vtkPoints()

    for i in range(n):
        multiplier = np.random.uniform(1, radius_max)

        point = voronoi.GetPoint(i)

        points.InsertNextPoint(point)
        cell_array.InsertNextCell(1)
        cell_array.InsertCellPoint(i)
        value = radius_array_data(i)
        radius_array.SetTuple(i, [value * multiplier])

    new_voronoi.SetPoints(points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)

    return new_voronoi


def read_command_line_surface(input_path=None, output_path=None):
    """
    Read arguments from commandline and return all values in a dictionary.
    If input_path and output_path are not None, then do not parse command line, but
    only return default values.
     Args:
        input_path (str): Input file path, positional argument with default None.
        output_path (str): Output file path, positional argument with default None.
    """
    description = " Remove or add noise within a specific relative or absolute bandwidth"
    parser = ArgumentParser(description=description, formatter_class=RawDescriptionHelpFormatter)

    # Add common arguments
    required = not (input_path is not None and output_path is not None)
    add_common_arguments(parser, required=required)

    # Set region of interest:
    parser.add_argument("-r", "--region-of-interest", type=str, default="full_model",
                        choices=["manual", "commandline", "first_line", "full_model"],
                        help="The method for defining the region to be changed. There are" +
                             " four options: 'manual', 'commandline', 'first_line', and" +
                             "'full'. In 'manual' the user will be provided with a" +
                             " visualization of the input surface, and asked to provide an" +
                             " end and start point of the region of interest." +
                             " If 'commandline' is provided, then '--region-points'" +
                             " is expected to be provided. In 'first_line' the section" +
                             " from the inlet to the first bifurcation will be chosen" +
                             " automatically. In 'full_model' the entire geometry will be chosen")
    parser.add_argument("--region-points", nargs="+", type=float, default=None, metavar="points",
                        help="If -r or --region-of-interest is 'commandline' then this" +
                             " argument have to be given. The method expects two points" +
                             " which defines the start and end of the region of interest. If" +
                             " 'method' is set to stenosis, then one point can be provided as well," +
                             " which is assumed to be the center of a new stenosis." +
                             " Example providing the points (1, 5, -1) and (2, -4, 3):" +
                             " --stenosis-points 1 5 -1 2 -4 3")

    parser.add_argument("-ns", "--noise", type=str2bool, default=False, metavar="Noise",
                        help="Turn on/off adding noise to the surface")

    parser.add_argument("-u", "--upper", type=float, default=1.0,
                        metavar="Upper smoothing factor",
                        help="Upper threshold for the smoothing")
    parser.add_argument("-fq", "--frequency", type=float, default=10, metavar="Frequency",
                        help="How frequent to add the noise. This number is drawn for each " +
                             " point along the centerline")
    parser.add_argument("-sd", "--frequency_deviation", type=float, default=2.5,
                        metavar="Standard deviation of frequency",
                        help="Standard deviation of distribution to draw frequency from.")
    parser.add_argument("-a", "--absolute", type=str2bool, default=False,
                        metavar="Absolute value", help="Turn on/off setting an absolute" +
                                                       " value for the smoothing criteria")
    parser.add_argument("-rm", "--radius-max", type=float, default=1.5,
                        metavar="Maximum radius", help="Draw a number from 1.0 and" +
                                                       " 'radius-max' and multiply with r to create the noise")

    # Output file argument
    if required:
        args = parser.parse_args()
    else:
        args = parser.parse_args(["-i" + input_path, "-o" + output_path])

    if args.region_of_interest == "commandline" and "region_points" is not None:
        raise ArgumentTypeError("When setting region of interest to commandline you have to" +
                                " provide the points through the argument '--regions-points'")

    # Do nothing warning
    if not args.smooth and not args.noise:
        print("WARNING: You have turned off both smoothing and adding noise. The algorithm" +
              " will therefore not do anything to the Voronoi diagram.")

    # Check smoothing_factor and upper values when not using absolute criteria
    if not args.absolute:
        if not 0 <= args.smooth_factor <= 1:
            raise ArgumentTypeError("When not using absolute values the 'smooth-factor'" + \
                                    "limit should be within [0, 1], not" + \
                                    " {}".format(args.smooth_factor()))
        if not 0 <= args.upper <= 1:
            raise ArgumentTypeError("When not using absolute values, the 'upper'" + \
                                    "limit should be within [0, 1], not" + \
                                    " {}".format(args.smooth_factor()))

    return dict(input_filepath=args.ifile, smooth=args.smooth, output_filepath=args.ofile,
                smooth_factor=args.smooth_factor, resampling_step=args.resampling_step,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point,
                poly_ball_size=args.poly_ball_size,
                region_of_interest=args.region_of_interest,
                region_points=args.region_points, upper=args.upper,
                frequency=args.frequency, noise=args.noise,
                frequency_deviation=args.frequency_deviation, absolute=args.absolute,
                radius_max=args.radius_max)


def main_surface():
    manipulate_surface(**read_command_line_surface())


if __name__ == "__main__":
    manipulate_surface(**read_command_line_surface())

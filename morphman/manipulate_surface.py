##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from argparse import ArgumentParser, RawDescriptionHelpFormatter

from morphman.common.argparse_common import *
from morphman.common.surface_operations import *
from morphman.common.vessel_reconstruction_tools import *


# Local import


def manipulate_surface(input_filepath, output_filepath, smooth, smooth_factor, no_smooth,
                       no_smooth_point, poly_ball_size, resampling_step,
                       region_of_interest, region_points, add_noise_lower_limit,
                       add_noise_upper_limit, frequency, frequency_deviation, noise,
                       absolute, radius_max, radius_min):
    """
    Controll the surface roughness by removing or adding points from the Voronoi diagram
    of the surface. A less general version of the smoothing algorithm was first presented
    in Ford et al. (2009).
     Args:
        input_filepath (str): Path to input surface.
        output_filepath (str): Path to output surface.
        smooth (bool): Determine if the voronoi diagram should be smoothed.
        smooth_factor (float): Smoothing factor used for voronoi diagram smoothing.
        no_smooth (bool): True of part of the model is not to be smoothed.
        no_smooth_point (ndarray): Point which is untouched by smoothing.
        poly_ball_size (list): Resolution of polyballs used to create surface.
        resampling_step (float): Resampling step used to resample centerlines.
        region_of_interest (str): Method for setting the region of interest ['manual' | 'commandline' | 'first_line' | 'full_model']
        region_points (list): If region_of_interest is 'commandline', this a flatten list of the start and endpoint
        add_noise_lower_limit (float): Upper bound for adding noise to the surface.
        add_noise_upper_limit (float): Upper bound for adding noise to the surface.
        frequency (float): Mean number of points added to the Voronoi diagram at each centerlinepoint.
        frequency_deviation (float): Standard deviation of frequency
        noise (bool): Turn on/off adding noise to the surface.
        absolute (bool): Absolute value for the smoothing criteria
        radius_max (float): Used to pick MISR multiplier to create noise on surface
    """
    # Filenames
    base_path = get_path_names(input_filepath)

    # Output filepaths
    centerlines_path = base_path + "_centerline.vtp"
    centerline_spline_path = base_path + "_centerline_region_of_interest.vtp"
    centerline_remaining_path = base_path + "_centerline_remaining.vtp"
    voronoi_new_path = base_path + "_voronoi_{}.vtp"
    if noise and smooth:
        voronoi_new_path = voronoi_new_path.format("smooth_noise")
    elif noise:
        voronoi_new_path = voronoi_new_path.format("noise")
    elif smooth:
        voronoi_new_path = voronoi_new_path.format("smooth")
    else:
        voronoi_new_path = voronoi_new_path.format("unchanged")

    # Clean and capp / uncapp surface
    surface, capped_surface = prepare_surface(base_path, input_filepath)

    # Get inlet and outlets
    inlet, outlets = get_inlet_and_outlet_centers(surface, base_path)
    centerlines, voronoi, pole_ids = compute_centerlines(inlet, outlets, centerlines_path,
                                                         capped_surface, resampling=resampling_step,
                                                         smooth=False, base_path=base_path)

    centerline_splined, centerline_remaining, centerline_diverging, region_points, _ = get_line_to_change(
        capped_surface, centerlines, region_of_interest, "bend", region_points, 0)
    write_polydata(centerline_splined, centerline_spline_path)

    # Split the Voronoi diagram
    if centerline_diverging is not None:
        centerline_remaining = vtk_merge_polydata([centerline_diverging] + centerline_remaining)
        write_polydata(centerline_remaining, centerline_remaining_path)

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

        voronoi_relevant = smooth_voronoi_diagram(voronoi_relevant, centerline_splined,
                                                  smooth_factor, no_smooth_cl, absolute)

    if noise:
        print("-- Add noise to the Voronoi diagram.")
        voronoi_relevant = add_noise_to_voronoi_diagram_new_points(surface,
                                                                   voronoi_relevant,
                                                                   centerline_splined,
                                                                   radius_max, radius_min,
                                                                   frequency,
                                                                   frequency_deviation,
                                                                   add_noise_lower_limit,
                                                                   add_noise_upper_limit,
                                                                   absolute)

    # Merge changed and unchanged Voronoi diagram
    voronoi = vtk_merge_polydata([voronoi_relevant, voronoi_rest])
    write_polydata(voronoi, voronoi_new_path)

    # Write a new surface from the new voronoi diagram
    print("-- Create new surface.")
    new_surface = create_new_surface(voronoi, poly_ball_size)

    print("-- Preparing surface for output.")
    new_surface = prepare_output_surface(new_surface, surface, centerlines,
                                         output_filepath, test_merge=False, changed=True,
                                         old_centerline=centerlines)

    print("-- Writing new surface to {}.".format(output_filepath))
    write_polydata(new_surface, output_filepath)


def add_noise_to_voronoi_diagram_new_points(surface, voronoi, centerline, radius_max,
                                            radius_min, frequency, frequency_deviation,
                                            lower, upper, absolute):
    """
    Add noise to Voronoi diagram by adjusting
    the MISR size by a factor in [1.0, radius_max],
    combined with a set frequency + deviation.

    Args:
        noise_method (string): Noise method which is applied to voronoi diagram
        voronoi (vtkPolyData): Voronoi Diagram to be smoothed
        centerline (vtkPolyData): Centerline(s) along relevant Voronoi diagram
        radius_max (float): Draw a number from 'radius-min' and 'radius-max'
                            and multiply with the distance from the new point
                            to the surface to set the radius of the new pointUsed
                            to pick MISR multiplier to create noise on surface.
        radius_min (float): Draw a number from 'radius-min' and 'radius-max'
                            and multiply with the distance from the new point
                            to the surface to set the radius of the new pointUsed
                            to pick MISR multiplier to create noise on surface
        frequency (float): Frequency at which noise is added to the voronoi diagram,
                           based on points along the centerline. Drawn from a normal
                           distribution with mean of 'frequency'.
        frequency_deviation (float): Standard deviation of frequency distribution.
        lower (float): The new location of a point in the voronoi diagram is a length
                       drawn between 'lower' and 'upper', either relative to the MISR
                       or in absolute value.
        upper (float): The new location of a point in the voronoi diagram is a length
                       drawn between 'lower' and 'upper', either relative to the MISR
                       or in absolute value.
        absolute (bool): Turn on/off an absolute threshold on the values instead of
                         relative with respect to the MISR.

    Returns:
        new_voronoi (vtkPolyData): Noisy Voronoi diagram
    """
    # Compute local coordinate system for each centerline point
    centerline = vmtk_compute_geometric_features(centerline, smooth=True)

    # Locator on the surface
    locator = get_vtk_point_locator(surface)

    # Boundaries
    boundaries = vtk_extract_feature_edges(surface)
    locator_boundaries = get_vtk_point_locator(boundaries)

    # Pointers to arrays
    frenet_tangent_data = centerline.GetPointData().GetArray("FrenetTangent").GetTuple1
    frenet_normal_data = centerline.GetPointData().GetArray("FrenetNormal").GetTuple1
    misr_data = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1

    # Number of new points
    n = centerline.GetNumberOfPoints()
    new_voronoi_points = np.rint(np.random.normal(frequency, frequency_deviation, n)).astype(int)
    new_voronoi_points[new_voronoi_points < 0] = 0
    number_of_new_points = np.sum(new_voronoi_points)

    # Holding the noise
    new_voronoi = vtk.vtkPolyData()
    cell_array = vtk.vtkCellArray()
    points = vtk.vtkPoints()
    radius_array = []

    # Add noise
    counter = 0
    for i in range(n):
        point = np.array(centerline.GetPoint(i))
        tangent = centerline.GetPointData().GetArray("FrenetTangent").GetTuple3(i)
        normal = centerline.GetPointData().GetArray("FrenetNormal").GetTuple3(i)
        misr = misr_data(i)

        for _ in range(new_voronoi_points[i]):
            # Define location of new point
            angle = np.random.uniform(0, 360)
            translation = vtk.vtkTransform()
            translation.RotateWXYZ(angle, tangent)
            tmp_normal = [0, 0, 0]
            translation.TransformNormal(normal, tmp_normal)
            if absolute:
                length = np.random.uniform(lower, upper)
            else:
                length = misr*np.random.uniform(lower, upper)
            new_point = point + length*np.array(tmp_normal)

            # Set new r
            surface_id = locator.FindClosestPoint(new_point)
            distance_to_surface = get_distance(surface.GetPoint(surface_id), new_point)
            new_r = distance_to_surface*np.random.uniform(radius_min, radius_max)

            # Check distance to the in- and outlets.
            boundaries_id = locator_boundaries.FindClosestPoint(new_point)
            distance_to_inlet = get_distance(boundaries.GetPoint(boundaries_id), new_point)
            if distance_to_inlet < new_r:
                continue

            # Add new point
            points.InsertNextPoint(new_point)
            cell_array.InsertNextCell(1)
            cell_array.InsertCellPoint(counter)
            radius_array.append(new_r)

            counter += 1

    radius_array = create_vtk_array(radius_array, radiusArrayName)
    new_voronoi.SetPoints(points)
    new_voronoi.SetVerts(cell_array)
    new_voronoi.GetPointData().AddArray(radius_array)

    new_voronoi = vtk_merge_polydata([new_voronoi, voronoi])

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

    # Variables for both smoothing surface and adding noise
    parser.add_argument("-a", "--absolute", type=str2bool, default=False,
                        metavar="Absolute value",
                        help="Turn on/off setting an absolute value for the smoothing" +
                             " criteria")

    # Variables for adding noise
    parser.add_argument("-ns", "--noise", type=str2bool, default=False, metavar="Noise",
                        help="Turn on/off adding noise to the surface")
    parser.add_argument("-u", "--add-noise-upper-limit", type=float, default=1.0,
                        metavar="Upper noise limit",
                        help="When adding noise the center of the added points are" +
                             " within in [lower, upper] * MISR from the center. A" +
                             " bandwidth on 0.9 and 1.0 would only provide rather" +
                             " 'high-frequent' noise on the surface.")
    parser.add_argument("-l", "--add_noise_lower_limit", type=float, default=0.85,
                        metavar="Lower noise limit.",
                        help="When adding noise the center of the added points are" +
                             " within in [lower, upper] * MISR from the center. A" +
                             " bandwidth on 0.9 and 1.0 would only provide rather" +
                             " 'high-frequent' noise on the surface.")
    parser.add_argument("-fq", "--frequency", type=float, default=5, metavar="Frequency",
                        help="How frequent to add the noise. This number is drawn for each " +
                             " point along the centerline")
    parser.add_argument("-sd", "--frequency-deviation", type=float, default=5,
                        metavar="Standard deviation of frequency",
                        help="Standard deviation of distribution to draw frequency from.")
    parser.add_argument("-rma", "--radius-max", type=float, default=1.5,
                        metavar="Maximum radius",
                        help="Draw a number from 'radius-min' and 'radius-max' and multiply with r" +
                             " to create the noise")
    parser.add_argument("-rmi", "--radius-min", type=float, default=1.0,
                        metavar="Maximum radius",
                        help="Draw a number from 'radius-min' and 'radius-max' and multiply with" +
                             " the distance from the new point to the surface to set the" +
                             " radius of the new point")


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
        if not 0 <= args.add_noise_lower_limit <= 1:
            raise ArgumentTypeError("When not using absolute values the 'smooth-factor'" +
                                    " limit should be within [0, 1], not" +
                                    " {}".format(args.smooth_factor()))
        if not 0 <= args.add_noise_lower_limit <= 1:
            raise ArgumentTypeError("When not using absolute values the 'add-noise-lower-limit'" +
                                    " limit should be within [0, 1], not" +
                                    " {}".format(args.add_noise_lower_limit()))
        if not 0 <= args.add_noise_upper_limit <= 1:
            raise ArgumentTypeError("When not using absolute values, the 'add-noise-lower-limit'" +
                                    " limit should be within [0, 1], not" +
                                    " {}".format(args.add_noise_upper_limit()))

    if args.frequency_deviation <= 0:
        raise ArgumentTypeError("The standard deviation has to be larger than zero." +
                                " Please provide a valid value, not {}"\
                                .format(args.frequency_deviation))

    return dict(input_filepath=args.ifile, smooth=args.smooth, output_filepath=args.ofile,
                smooth_factor=args.smooth_factor, resampling_step=args.resampling_step,
                no_smooth=args.no_smooth, no_smooth_point=args.no_smooth_point,
                poly_ball_size=args.poly_ball_size,
                region_of_interest=args.region_of_interest,
                region_points=args.region_points,
                add_noise_upper_limit=args.add_noise_upper_limit,
                add_noise_lower_limit=args.add_noise_lower_limit,
                frequency=args.frequency, noise=args.noise,
                frequency_deviation=args.frequency_deviation, absolute=args.absolute,
                radius_max=args.radius_max, radous_min=args.radius_min)


def main_surface():
    manipulate_surface(**read_command_line_surface())


if __name__ == "__main__":
    manipulate_surface(**read_command_line_surface())

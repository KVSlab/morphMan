##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.


# This demo is the equivalent of executing the following in a terminal:
#  $ morphman-surface --ifile C0005/surface/model.vtp --ofile C0005/surface/model_surface_noise.vtp \
#       --poly-ball-size 250 250 250
#
#  $ morphman-surface --ifile C0005/surface/model.vtp --ofile C0005/surface/model_surface_smooth.vtp \
#       --noise True --smooth False --fq 5 -sd 0.1 -u 0.95 -l 0.80 -rma 1.1 -rmi 1.35 \
#       --poly-ball-size 250 250 250
#
# and a more detailed explenation could be found here:
# https://morphman.readthedocs.io/en/latest/manipulate_surface.html#tutorial-manipulate-surface-roughness

from morphman import manipulate_surface, main_surface, read_command_line_surface
from os import path
from get_test_data import download_case

# Set absolute path to the demo folder
absolute_path = path.dirname(path.abspath(__file__))

# Set input and output paths
case = "C0005"
input_filepath = path.join(absolute_path, case, "surface", "model.vtp")
output_filepath = path.join(absolute_path, case, "surface", "surface_smooth.vtp")

# Download case from the Aneurisk web for this demo
if not path.exists(path.join(input_filepath)):
    download_case(case)

# Get default values for manipulate_surface
default_values = read_command_line_surface(input_filepath, output_filepath)

# Problem specific parameters
default_values["no_smooth"] = True
default_values["no_smooth_point"] = [54.7310791015625, 46.95508575439453, 41.97511291503906]

# Parameters for reconstructing the surface
default_values["poly_ball_size"] = [250, 250, 250]

# Run manipulation - smooth the entire geometry except aneurysm dome
manipulate_surface(**default_values)

### Add noise to the surface
# Output file path
default_values["output_filepath"] = path.join(absolute_path, case, "surface", "surface_noise.vtp")

# Method specific paramters
default_values["smooth"] = False
default_values["noise"] = True
default_values["frequency"] = 0
default_values["frequency_deviation"] = 1
default_values["add_noise_lower_limit"] = 0.8
default_values["add_noise_upper_limit"] = 0.90
default_values["radius_min"] = 1.1
default_values["radius_max"] = 1.5

manipulate_surface(**default_values)

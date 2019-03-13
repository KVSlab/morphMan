##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.


# This demo is the equivalent of executing the following in a terminal:
#  $ morphman-curvature --ifile C0005/surface/model.vtp --ofile C0005/surface/model_curvature_decreased.vtp \
#       --smooth-line True --iterations 100 --smooth-factor-line 1.8 --region-of-interest first_line \
#       --poly-ball-size 250 250 250
#
#  $ morphman-curvature --ifile C0005/surface/model.vtp --ofile C0005/surface/model_curvature_increased.vtp \
#       --smooth-line False --iterations 100 --smooth-factor-line 1.8 --region-of-interest first_line \
#       --poly-ball-size 250 250 250
#
# and a more detailed explenation could be found here:
# https://morphman.readthedocs.io/en/latest/manipulate_curvature.html#tutorial-manipulate-curvature

from morphman import manipulate_curvature, main_curvature, read_command_line_curvature
from os import path
from get_test_data import download_case

# Set absolute path to the demo folder
absolute_path = path.dirname(path.abspath(__file__))

# Set input and output paths
case = "C0005"
input_filepath = path.join(absolute_path, case, "surface", "model.vtp")
output_filepath = path.join(absolute_path, case, "surface", "curvature_plus.vtp")

# Download case from the Aneurisk web for this demo
if not path.exists(path.join(input_filepath)):
    download_case(case)

# Get default values for manipulate_curvature
default_values = read_command_line_curvature(input_filepath, output_filepath)

# Set region of interest
default_values["region_of_interest"] = "first_line"

# Method spesific parameters - reduce the curvature
default_values["smooth_line"] = True
default_values["iterations"] = 100
default_values["smooth_factor_line"] = 1.8

# Parameters for reconstructing the surface
default_values["poly_ball_size"] = [250, 250, 250]

# Run manipulation
manipulate_curvature(**default_values)


### Increase the curvature
# Method specific paramters
default_values["smooth_line"] = False

# Set new output path
output_filepath = path.join(absolute_path, case, "surface", "no_noch.vtp")

manipulate_curvature(**default_values)

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.


# This demo is the equivalent of executing the following in a terminal:
#  $ morphman-bend --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_vertical_plus.vtp \
#       --alpha 0.4 --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 \
#       --poly-ball-size 250 250 250
#  $ morphman-bend --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_vertical_minus.vtp \
#       --alpha -0.4 --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 \
#       --poly-ball-size 250 250 250
#
#  $ morphman-bend --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_horizontal_plus.vtp \
#       --beta 0.4 --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 \
#       --poly-ball-size 250 250 250
#  $ morphman-bend --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_horizontal_minus.vtp \
#       --beta -0.4 --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 \
#       --poly-ball-size 250 250 250
#
#  $ morphman-bend --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_plus.vtp \
#       --alpha 0.4 --beta 0.4 --region-of-interest commandline \
#       --region-points 49.8 49.7 36.6 53.1 41.8 38.3 --poly-ball-size 250 250 250
#  $ morphman-bend --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_minus.vtp \
#       --alpha -0.4 --beta -0.4 --region-of-interest commandline \
#       --region-points 49.8 49.7 36.6 53.1 41.8 38.3 --poly-ball-size 250 250 250
#
# and a more detailed explenation could be found here:
# https://morphman.readthedocs.io/en/latest/manipulate_bend.html#tutorial-manipulate-bend

from morphman import manipulate_bend, main_bend, read_command_line_bend
from os import path
from get_test_data import download_case

# Set absolute path to the demo folder
absolute_path = path.dirname(path.abspath(__file__))

# Set input and output paths
case = "C0005"
input_filepath = path.join(absolute_path, case, "surface", "model.vtp")
output_filepath = path.join(absolute_path, case, "surface", "%s.vtp")

# Download case from the Aneurisk web for this demo
if not path.exists(path.join(input_filepath)):
    download_case(case)

# Get default values for manipulate_bend
default_values = read_command_line_bend(input_filepath, output_filepath)

# Set region of interest
default_values["region_of_interest"] = "commandline"
default_values["region_points"] = [49.8, 49.7, 36.6, 53.1, 41.8, 38.3]

# Method spesific parameters
beta = [0, 0, 0.4, -0.4, 0.4, -0.4]
alpha = [0.4, -0.4, 0, 0, 0.4, -0.4]
names = ["bend_vertical_plus", "bend_vertical_minus", "bend_horizonal_plus",
         "bend_horizontal_plus", "bend_plus", "bend_minus"]

# Parameters for reconstructing the surface
default_values["poly_ball_size"] = [250, 250, 250]

# Run manipulation
for b, a, n in zip(beta, alpha, names):
    default_values["beta"] = b
    default_values["alpha"] = a
    default_values["output_filepath"] = n + ".vtp"
    manipulate_bend(**default_values)

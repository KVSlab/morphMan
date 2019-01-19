##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.


# This demo is the equivalent of executing the following in a terminal:
#  $ morphman-area --ifile C0002/surface/model.vtp --ofile C0002/surface/stenosis.vtp \
#       --method stenosis --size 4 --percentage 50 --region-of-interest commandline \
#       --region-points 28.7 18.4 39.5 --poly-ball-size 250 250 250
#  $ morphman-area --ifile C0002/surface/stenosis.vtp --ofile C0002/surface/stenosis_removed.vtp \
#       --method stenosis --region-of-interest commandline \
#       --region-points 30.1 18.5 34.6 27.1 12.7 38.2 --poly-ball-size 250 250 250
#
# and a more detailed explenation could be found here:
# https://morphman.readthedocs.io/en/latest/manipulate_area.html#create-remove-a-stenosis

from morphman import manipulate_area, main_area, read_command_line_area
from os import path
from get_test_data import download_case

# Set absolute path to the demo folder
absolute_path = path.dirname(path.abspath(__file__))

# Set input and output paths
case = "C0002"
input_filepath = path.join(absolute_path, case, "surface", "model.vtp")
output_filepath = path.join(absolute_path, case, "surface", "stenosis.vtp")

# Download case from the Aneurisk web for this demo
if not path.exists(path.join(input_filepath)):
    download_case(case)

# Get default values for manipulate_area
default_values = read_command_line_area(input_filepath, output_filepath)

# Set region of interest
default_values["region_of_interest"] = "commandline"
default_values["region_points"] = [28.7, 18.4, 39.5]

# Method for changing the area
default_values["method"] = "stenosis"

# Method spesific parameters - create a stenosis
default_values["size"] = 4.0
default_values["percentage"] = 50

# Parameters for reconstructing the surface
default_values["poly_ball_size"] = [250, 250, 250]

# Run manipulation
manipulate_area(**default_values)


### Decrease cross-section area ratio
# Method spesific parameters
default_values["region_points"] = [30.1, 18.5, 34.6, 27.1, 12.7, 38.2]

# Set new output path
output_filepath = path.join(absolute_path, case, "surface", "stenosis_removed.vtp")

# Run manipulation
manipulate_area(**default_values)

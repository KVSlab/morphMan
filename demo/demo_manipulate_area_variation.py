##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.


# This demo is the equivalent of executing the following in a terminal:
#  $ morphman-area --ifile C0002/surface/model.vtp --ofile C0002/surface/increased_variation.vtp \
#        --method variation --ratio 4.0 --region-of-interest first_line --poly-ball-size 250 250 250
#  $ morphman-area --ifile C0002/surface/model.vtp --ofile C0002/surface/decreased_variation.vtp \
#        --method variation --ratio 0.5 --region-of-interest first_line --poly-ball-size 250 250 250
#
# and a more detailed explenation could be found here:
# https://morphman.readthedocs.io/en/latest/manipulate_area.html#area-variations

from morphman import manipulate_area, main_area, read_command_line_area
from os import path
from get_test_data import download_case

# Set absolute path to the demo folder
absolute_path = path.dirname(path.abspath(__file__))

# Set input and output paths
case = "C0002"
input_filepath = path.join(absolute_path, case, "surface", "model.vtp")
output_filepath = path.join(absolute_path, case, "surface", "increased_variation.vtp")

# Download case from the Aneurisk web for this demo
if not path.exists(path.join(input_filepath)):
    download_case(case)

# Get default values for manipulate_area
default_values = read_command_line_area(input_filepath, output_filepath)

# Set region of interest
default_values["region_of_interest"] = "first_line"

# Method for changing the area
default_values["method"] = "variation"

# Method spesific parameters - increase cross-section area ratio
default_values["ratio"] = 4.0

# Parameters for reconstructing the surface
default_values["poly_ball_size"] = [250, 250, 250]

# Run manipulation
manipulate_area(**default_values)


### Decrease cross-section area ratio
# Method spesific parameters
default_values["ratio"] = 1.0

# Set new output path
output_filepath = path.join(absolute_path, case, "surface", "decreased_variation.vtp")

# Run manipulation
manipulate_area(**default_values)

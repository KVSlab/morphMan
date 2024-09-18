##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.
## This software is distributed WITHOUT ANY WARRANTY; without even
## the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE.  See the above copyright notices for more information.


# This demo is the equivalent of executing the following in a terminal:
#  $ morphman-bifurcation --ifile C0005/surface/model.vtp --ofile C0005/surface/rotate_plus.vtp \
#       --angle 20 --region-of-interest commandline --region-points 43.2 70.5 26.4 84.4 60.6 50.6 \
#       --poly-ball-size 250 250 250
#  $ morphman-bifurcation --ifile C0005/surface/model.vtp --ofile C0005/surface/rotate_minus.vtp \
#       --angle -20 --region-of-interest commandline --region-points 43.2 70.5 26.4 84.4 60.6 50.6 \
#        --poly-ball-size 250 250 250
#
#  $ morphman-bifurcation --ifile C0005/surface/model.vtp --ofile C0005/surface/rotate_no_notch.vtp \
#       --angle -20 --bif True --lower True --region-of-interest commandline \
#       --region-points 43.2 70.5 26.4 84.4 60.6 50.6 --poly-ball-size 250 250 250
#
#  $ morphman-bifurcation --ifile C0066/surface/model.vtp --ofile C0066/surface/removed_aneurysm.vtp \
#       --keep-fixed-1 True --keep-fixed-2 True --bif True --lower True --angle 0 \
#       --region-of-interest commandline --region-points 31.37 60.65 25.21 67.81 43.08 41.24 \
#       --poly-ball-size 250 250 250
#
# and a more detailed explenation could be found here:
# https://morphman.readthedocs.io/en/latest/manipulate_bifurcation.html#tutorial-manipulate-bifurcation

from os import path

from get_test_data import download_case

from morphman.manipulate_bifurcation import (
    manipulate_bifurcation,
    read_command_line_bifurcation,
)

# Set absolute path to the demo folder
absolute_path = path.dirname(path.abspath(__file__))

# Set input and output paths
case = "C0005"
input_filepath = path.join(absolute_path, case, "surface", "model.vtp")
output_filepath = path.join(absolute_path, case, "surface", "rotate_plus.vtp")

# Download case from the Aneurisk web for this demo
if not path.exists(path.join(input_filepath)):
    download_case(case)

# Get default values for manipulate_bifurcation
default_values = read_command_line_bifurcation(input_filepath, output_filepath)

# Set region of interest
default_values["region_of_interest"] = "commandline"
default_values["region_points"] = [43.2, 70.5, 26.4, 84.4, 60.6, 50.6]

# Method spesific parameters - rotate the daughter branches together
default_values["angle"] = 20

# Parameters for reconstructing the surface
default_values["poly_ball_size"] = [250, 250, 250]

# Run manipulation
manipulate_bifurcation(**default_values)

### Negative angle
# Method spesific parameters - rotate the daughter branches appart
default_values["angle"] = -20

# Set new output path
output_filepath = path.join(absolute_path, case, "surface", "rotate_minus.vtp")

# Run manipulation
manipulate_bifurcation(**default_values)

### No-noch
# Method spesific parameters - improved rebuilding of the bifurcation
default_values["bif"] = True
default_values["lower"] = True

# Set new output path
output_filepath = path.join(absolute_path, case, "surface", "no_noch.vtp")

# Run manipulation
manipulate_bifurcation(**default_values)

### Remove an aneurysm
case = "C0066"
default_values["input_filepath"] = path.join(
    absolute_path, case, "surface", "model.vtp"
)
default_values["output_filepath"] = path.join(
    absolute_path, case, "surface", "rotate_plus.vtp"
)

# Download case from the Aneurisk web for this demo
if not path.exists(path.join(default_values["input_filepath"])):
    download_case(case)

# Set region of interest
default_values["region_of_interest"] = "commandline"
default_values["region_points"] = [31.37, 60.65, 25.21, 67.81, 43.08, 41.24]

# Method spesific parameters
default_values["keep_fixed_1"] = True
default_values["keep_fixed_2"] = True
default_values["bif"] = True
default_values["lower"] = True
default_values["angle"] = 0

# Run manipulation
manipulate_bifurcation(**default_values)

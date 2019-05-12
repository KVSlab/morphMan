##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.


# This demo is the equivalent of executing the following in a terminal:
#   $ morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/moved_branch.vtp --branch-number 1
#       --branch-location 21.7 18.1 25.9  --translation-method commandline --poly-ball-size 250 250 250
#
#   $ morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/moved_and_rotated_branch.vtp
#       --azimuth-angle 180 --branch-number 1 --branch-location 21.7 18.1 25.9
#       --translation-method commandline  --poly-ball-size 250 250 250
#
#   $ morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/azimuth_rotated_branch.vtp
#       --azimuth-angle 120 --branch-number 1 --translation-method no_translation  --poly-ball-size 250 250 250
#
#   $ morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/polar_rotated_branch.vtp
#       --polar-angle 20 --branch-number 1 --translation-method no_translation  --poly-ball-size 250 250 250
#
#   $ morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/polar_rotated_branch.vtp
#       --polar-angle -20 --branch-number 1 --translation-method no_translation  --poly-ball-size 250 250 250
#
#   $ morphman-branch --ifile C0002/surface/model.vtp --ofile C0002/surface/removed_branch.vtp
#       --remove-branch True --branch-number 4 --poly-ball-size 250 250 250
#
# and a more detailed explenation could be found here:
# https://morphman.readthedocs.io/en/latest/manipulate_branch.html

from os import path
from get_test_data import download_case
from morphman import manipulate_branch, main_branch, read_command_line_branch

# Set absolute path to the demo folder
absolute_path = path.dirname(path.abspath(__file__))

# Set input and output paths
case = "C0002"
input_filepath = path.join(absolute_path, case, "surface", "model.vtp")
output_filepath = path.join(absolute_path, case, "surface", "moved_branch.vtp")

# Download case from the Aneurisk web for this demo
if not path.exists(path.join(input_filepath)):
    download_case(case)

# Get default values for manipulate_branch
default_values = read_command_line_branch(input_filepath, output_filepath)

### Translation of branch
# Set region of interest
default_values["branch_to_manipulate_number"] = 1
default_values["branch_location"] = [21.7, 18.1, 25.9]
default_values["translation_method"] = "commandline"

# Parameters for reconstructing the surface
default_values["poly_ball_size"] = [250, 250, 250]

# Set new output path
default_values["output_filepath"] = output_filepath

# Run manipulation
manipulate_branch(**default_values)

### Translation and rotation of branch
# Method spesific parameters - rotation and translation of branch
default_values["azimuth_angle"] = 180

# Set new output path
default_values["output_filepath"] = path.join(absolute_path, case, "surface", "moved_and_rotated_branch.vtp")

# Run manipulation
manipulate_branch(**default_values)


### Rotation around branch base
# Method spesific parameters - rotation only
default_values["azimuth_angle"] = 120
default_values["translation_method"] = "no_translation"

# Set new output path
default_values["output_filepath"] = path.join(absolute_path, case, "surface", "azimuth_rotated_branch.vtp")

# Run manipulation
manipulate_branch(**default_values)

### Rotation around surface tangent
# Method spesific parameters - rotation only
default_values["azimuth_angle"] = 0
default_values["polar_angle"] = 20

# Set new output path
default_values["output_filepath"] = path.join(absolute_path, case, "surface", "polar_rotated_branch_up.vtp")

### Rotation around surface tangent
# Method spesific parameters - rotation only
default_values["polar_angle"] = -20

# Set new output path
default_values["output_filepath"] = path.join(absolute_path, case, "surface", "polar_rotated_branch_down.vtp")

# Run manipulation
manipulate_branch(**default_values)

### Clamp branch when translating
# Reset angle parameters
default_values["polar_angle"] = 0
default_values["azimuth_angle"] = 0

# Method spesific parameters
default_values["branch_location"] = [27.75, 29.75, 24.61]
default_values["clamp_branch"] = True

# Set new output path
default_values["output_filepath"] = path.join(absolute_path, case, "surface", "clamped_branch.vtp")

# Run manipulation
manipulate_branch(**default_values)


### Branch removal
# Method spesific parameters
default_values["branch_to_manipulate_number"] = 4
default_values["remove_branch"] = True

# Set new output path
default_values["output_filepath"] = path.join(absolute_path, case, "surface", "removed_branch.vtp")

# Run manipulation
manipulate_branch(**default_values)

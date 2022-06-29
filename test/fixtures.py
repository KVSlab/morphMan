##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from os import system, path, makedirs
from sys import platform

import pytest

from morphman.common import get_inlet_and_outlet_centers, get_path_names, compute_centerlines, prepare_surface


def download_testdata(test_path, output_file):
    if platform == "darwin":
        system("curl -L {} --output {}".format(test_path, output_file))
    elif platform == "linux" or platform == "linux2":
        system("wget -O {} {}".format(output_file, test_path))
    elif platform == "win32":
        system("bitsadmin /transfer download_model /download /priority high {} {}".format(test_path, output_file))


@pytest.fixture(scope="module")
def surface_paths():
    abs_path = path.dirname(path.abspath(__file__))

    # Path to test data
    # TODO: Replace with official Aneurisk database when available
    test_path = "https://github.com/hkjeldsberg/AneuriskDatabase/raw/master/models/C0001/surface/model.vtp"
    output_file = path.join(abs_path, "C0001", "surface", "model.vtp")

    # Create test data folders
    if not path.exists(path.join(abs_path, "C0001", "surface")):
        makedirs(path.join(abs_path, "C0001", "surface"))

    # Download test data if necessary
    if not path.exists(path.join(abs_path, "C0001", "surface", "model.vtp")):
        try:
            download_testdata(test_path, output_file)
        except Exception:
            raise Exception("Problem downloading the testdata, please do it manually from "
                            + test_path + " and extract the compressed tarball in the test folder")

    # Create centerline
    input_filepath = path.join(abs_path, "C0001", "surface", "model.vtp")
    base_path = get_path_names(input_filepath)
    centerlines_path = base_path + "_centerline.vtp"
    if not path.exists(centerlines_path):
        surface, capped_surface = prepare_surface(base_path, input_filepath)
        inlet, outlets = get_inlet_and_outlet_centers(surface, base_path)
        compute_centerlines(inlet, outlets, centerlines_path, capped_surface, resampling=0.1, smooth=False,
                            base_path=base_path)

    return [input_filepath, input_filepath.replace(".vtp", "_output.vtp")]

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from os import system, path
from sys import platform
from morphman.common import get_centers, get_path_names, prepare_surface, compute_centerlines

import pytest


def download_testdata(test_path, outputfile):
    if platform == "darwin":
        system("curl {} --output {}".format(test_path, outputfile))
        system("tar -zxvf {}".format(outputfile))
        system("rm {}".format(outputfile))
    elif platform == "linux" or platform == "linux2":
        system("wget {}".format(test_path))
        system("tar -zxvf {}".format(outputfile))
        system("rm {}".format(outputfile))
    elif platform == "win32":
        system("bitsadmin /transfer download_model /download /priority high {} {}".format(test_path, outputfile))
        system("tar -zxvf {}".format(outputfile))
        system("del /f {}".format(outputfile))


@pytest.fixture(scope="module")
def surface_paths():
    abs_path = path.dirname(path.abspath(__file__))

    # Path to test data
    test_path = "http://ecm2.mathcs.emory.edu/aneuriskdata/download/C0001/C0001_models.tar.gz"
    outputfile = path.join(abs_path, "C0001_models.tar.gz")

    # Download test data if necessary
    if not path.exists(path.join(abs_path, "C0001")):
        try:
            download_testdata(test_path, outputfile)
        except Exception:
            raise Exception("Problem downloading the testdata, please do it manually from " \
                            + test_path + " and extract the compressed tarball in the" \
                            + " test folder")

    # Create centerline
    input_filepath = path.join(abs_path, "C0001", "surface", "model.vtp")
    base_path = get_path_names(input_filepath)
    centerlines_path = base_path + "_centerline.vtp"
    if not path.exists(centerlines_path):
        surface, capped_surface = prepare_surface(base_path, input_filepath)
        inlet, outlets = get_centers(surface, base_path)
        centerlines, _, _ = compute_centerlines(inlet, outlets, centerlines_path,
                                                capped_surface, resampling=0.1,
                                                smooth=False, base_path=base_path)

    return [input_filepath, input_filepath.replace(".vtp", "_output.vtp")]

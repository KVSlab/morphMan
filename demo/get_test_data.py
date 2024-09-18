##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.
## This software is distributed WITHOUT ANY WARRANTY; without even
## the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
## PURPOSE.  See the above copyright notices for more information.

from os import makedirs, path, system
from sys import platform


def download_case(case):
    abs_path = path.dirname(path.abspath(__file__))
    output_file = path.join(abs_path, case, "surface", "model.vtp")
    url = "https://github.com/hkjeldsberg/AneuriskDatabase/raw/master/models/{}/surface/model.vtp".format(
        case
    )

    # Create test data folders
    if not path.exists(path.join(abs_path, case, "surface")):
        makedirs(path.join(abs_path, case, "surface"))

    try:
        if platform == "darwin":
            system("curl -L {} --output {}".format(url, output_file))
        elif platform == "linux" or platform == "linux2":
            system("wget -O {} {}".format(output_file, url))
        elif platform == "win32":
            system(
                "bitsadmin /transfer download_model /download /priority high {} {}".format(
                    url, output_file
                )
            )

    except Exception:
        raise RuntimeError(
            "Problem downloading the testdata, please do it manually from "
            + url
            + " and extract the compressed tarball in the"
            + " test folder"
        )


if __name__ == "__main__":
    download_case("C0001")
    download_case("C0002")
    download_case("C0003")
    download_case("C0004")
    download_case("C0005")

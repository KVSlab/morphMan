##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even 
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
##      PURPOSE.  See the above copyright notices for more information.

from os import system, path
from sys import platform


def download_case(case):
    abs_path = path.dirname(path.abspath(__file__))
    output_file = path.join(abs_path, "{}_models.tar.gz".format(case))
    adress = "http://ecm2.mathcs.emory.edu/aneuriskdata/download/{}/{}_models.tar.gz".format(case, case)

    try:
        if platform == "darwin":
            system("curl {} --output {}".format(adress, output_file))
            system("tar -zxvf {}".format(output_file))
            system("rm {}".format(output_file))

        elif platform == "linux" or platform == "linux2":
            system("wget {}".format(adress))
            system("tar -zxvf {}".format(output_file))
            system("rm {}".format(output_file))

        elif platform == "win32":
            system("bitsadmin /transfer download_model /download /priority high {} {}".format(adress, output_file))
            system("tar -zxvf {}".format(output_file))
            system("del /f {}".format(output_file))

    except:
        raise RuntimeError("Problem downloading the testdata, please do it manually from "
                           + adress + " and extract the compressed tarball in the"
                           + " test folder")


if __name__ == "__main__":
    download_case("C0001")
    download_case("C0002")
    download_case("C0003")
    download_case("C0004")
    download_case("C0005")

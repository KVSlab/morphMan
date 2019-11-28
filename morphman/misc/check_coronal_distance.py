##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter

# Local import
from morphman.common import *
import matplotlib.pyplot as plt


def iterate_models():
    dir = os.path.dirname(os.path.abspath(__file__))
    cases = os.listdir(dir + "/cases")
    resampling_step = 0.1
    show_plot = True
    counter = 1
    for case in cases:
        input_filepath = "./%s/surface/model.vtp"
        model = read_polydata(input_filepath)
        base_path = get_path_names(input_filepath)

        # Extract carotid siphon
        line = extract_ica_centerline(base_path, input_filepath, resampling_step)
        length = get_curvilinear_coordinate(line)

        x = np.zeros(length.shape[0])
        y = np.zeros(length.shape[0])
        z = np.zeros(length.shape[0])
        for i in range(z.shape[0]):
            x[i] = line.GetPoints().GetPoint(i)[0]
            y[i] = line.GetPoints().GetPoint(i)[1]
            z[i] = line.GetPoints().GetPoint(i)[2]

        minus_length = -1 * length
        if show_plot:
            plt.plot(length, z)
            plt.show()

        start = z[0]
        end = z[-1]

        if start > end:
            print("Coronal max is found at inlet, n=%i" % counter)
            counter += 1


iterate_models()

##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.
from IPython import embed

from morphman.common import *

import json, os


def test_get_path_names():
    filepath = "surface/model.vtp"
    base_path = get_path_names(filepath)

    assert base_path == "surface/model"


def test_get_distance():
    p1 = np.array([-1, -1, -1])
    p2 = np.array([1, 1, 1])

    distance = get_distance(p1, p2)
    true_distance = 2 * np.sqrt(3)

    assert distance == true_distance


def test_gram_schmidth():
    A = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]])
    E = gram_schmidt(A)

    e1 = E[0]
    e2 = E[1]
    e3 = E[2]

    products = np.abs([e1.dot(e2), e2.dot(e3), e1.dot(e3)])
    tol = 1e-15

    assert sum(products) < tol


def test_get_parameters():
    # Add mock parameters

    params = {'key': 'value'}

    with open('test_info.json', 'w') as outfile:
        json.dump(params, outfile)

    # Test
    data = get_parameters("test")

    assert data['key'] == 'value'

    # Remove mock data
    os.remove("test_info.json")


def test_write_parameters():
    params = {'key': 'value'}
    write_parameters(params, "test")

    assert os.path.isfile("test_info.json")

    os.remove("test_info.json")


def test_numpy_to_polydata():
    v = np.array([[0, 0, 0], [1, 1, 1], [2, 2, 2]])

    v_as_vtk_array = convert_numpy_data_to_polydata(v, ["zero", "one", "two"])

    assert v_as_vtk_array.GetPoint(0) == (0, 0, 0)
    assert v_as_vtk_array.GetPoint(1) == (1, 1, 1)
    assert v_as_vtk_array.GetPoint(2) == (2, 2, 2)

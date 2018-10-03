import pytest
from manipulate_area import *

relative_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(relative_path, '..', 'src'))


@pytest.fixture(scope='module')
def init_data():
    # Global parameters
    input_filepath = "./testdata/P0134/surface/model.vtp"
    output_filepath = "./testdata/P0134/surface/model_output_area_stenosis.vtp"
    smooth_factor = 0.25
    poly_ball_size = [120, 120, 120]
    smooth = True
    no_smooth = False
    no_smooth_point = None
    method = "area"
    stenosis_length = 2.0
    ratio = None
    beta = 0.8
    percentage = 20

    return dict(input_filepath=input_filepath, method=method, smooth=smooth,
                smooth_factor=smooth_factor, percentage=percentage,
                ratio=ratio, stenosis_length=stenosis_length,
                beta=beta, output_filepath=output_filepath,
                poly_ball_size=poly_ball_size, no_smooth=no_smooth,
                no_smooth_point=no_smooth_point)


def test_create_stenosis(init_data):
    # Get region points
    base_path = get_path_names(init_data['input_filepath'])
    centerline = extract_single_line(read_polydata(base_path + "_centerline.vtp"), 0)
    n = centerline.GetNumberOfPoints()
    region_point = list(centerline.GetPoint(int(n * 0.4)))
    init_data['region_of_interest'] = "commandline"
    init_data['region_points'] = region_point

    # Perform area manipulation
    area_variations(**init_data)

    # Import old area and splined centerline for region of interest
    base_path = get_path_names(init_data['input_filepath'])
    centerline_spline_path = base_path + "_centerline_spline.vtp"
    centerline_area_spline_path = base_path + "_centerline_area_spline.vtp"
    new_surface_path = base_path + "_output_area_pure_variation.vtp"
    surface = read_polydata(new_surface_path)
    centerline_area = read_polydata(centerline_area_spline_path)
    centerline_spline = read_polydata(centerline_spline_path)
    new_centerline_area, _ = vmtk_compute_centerline_sections(surface,
                                                              centerline_spline)
    old_area = get_array("CenterlineSectionArea", centerline_area)
    new_area = get_array("CenterlineSectionArea", new_centerline_area)
    area_decreased = []

    # Ignore ends as these are adjusted less / area computation may contain errors
    start = int(old_area.size * 0.05)
    stop = int(old_area.size * 0.95)
    for i in range(start, stop):
        area_decreased.append(old_area[i] > new_area[i])

    assert area_decreased.count(True) == 0

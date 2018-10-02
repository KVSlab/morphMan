import pytest
from manipulate_area import *

relative_path = path.dirname(path.abspath(__file__))
sys.path.insert(0, path.join(relative_path, '..', 'src'))


@pytest.fixture(scope='module')
def init_data():
    # Global parameters
    input_filepath = "./testdata/P0134/surface/model.vtp"
    output_filepath = "./testdata/P0134/surface/model_output_area_variation.vtp"
    smooth_factor = 0.25
    poly_ball_size = [120, 120, 120]
    smooth = True
    no_smooth = False
    no_smooth_point = None
    method = "variation"
    region_of_interest = "first_line"
    region_points = None
    stenosis_length = 2.0
    percentage = 50
    ratio = None

    return dict(input_filepath=input_filepath, method=method, smooth=smooth,
                smooth_factor=smooth_factor,
                region_of_interest=region_of_interest,
                region_points=region_points, ratio=ratio,
                stenosis_length=stenosis_length,
                percentage=percentage, output_filepath=output_filepath,
                poly_ball_size=poly_ball_size, no_smooth=no_smooth,
                no_smooth_point=no_smooth_point)


def test_increased_ratio(init_data):
    # Increase R (ratio)
    init_data['beta'] = 0.8
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
    mean_area = np.mean(old_area)
    variations = []
    # FIXME:  What is the expected outcome?
    for i in range(old_area.size):
        if old_area[i] <= mean_area:
            variations.append(old_area[i] >= new_area[i])
        else:
            variations.append(old_area[i] < new_area[i])
    assert variations.count(False) == 0

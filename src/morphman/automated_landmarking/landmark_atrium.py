## Adhoc landmarking script of the left atrium
## method is a generalized version of the method found in:
## Tobon-Gomez, Catalina, et al. "Benchmark for algorithms segmenting
## the left atrium from 3D CT and MRI datasets." IEEE transactions on medical
## imaging 34.7 (2015): 1460-1473.

## Test data can be aquired from the Left Atrium Segmentation Challenge
## http://www.cardiacatlas.org/challenges/left-atrium-segmentation-challenge/

## Writen by Aslak W. Bergersen, 2019
## Modified by Henrik A. Kjeldsberg, 2022

from morphman.common import *


def extract_LA_and_LAA():
    """Algorithm for detecting the left atrial appendage and isolate it from the atrium lumen
     based on the cross-sectional area along enterlines.

    Args:
        file_path (str): Path to the surface for landmarking
        input_path (str): Path to store the landmarked surface

    Output:
        surface (vtkPolyData): A landmarked surface
    """

    # File paths
    print("--- Load model file\n")
    save_path = "/Users/hkjeldsberg/Simula/morphMan/models"
    input_path = path.join(save_path, "model.vtp")
    clipped_model = input_path.replace(".vtp", "_la_and_laa.vtp")

    new_cl = input_path.replace(".vtp", "_centerline.vtp")
    surface = read_polydata(input_path)

    capped_surface = vmtk_cap_polydata(surface)

    # Centers
    inlet, outlets = compute_centers(surface)
    p_outlet = np.array(inlet)

    # Centerline
    capped_surface = vtk_clean_polydata(capped_surface)
    la_centerlines, _, _ = compute_centerlines(inlet, outlets, new_cl, capped_surface, resampling=0.1, smooth=True)

    # Clip PVs
    N_PVs = 4
    for i in range(N_PVs):
        print(f"--- Clipping PV ({i + 1})")
        la_centerline_i = extract_single_line(la_centerlines, i)
        start = int(la_centerline_i.GetNumberOfPoints() * 0.5)
        stop = int(la_centerline_i.GetNumberOfPoints() * 0.95)

        line = extract_single_line(la_centerline_i, 0, start_id=start, end_id=stop)
        la_l = get_curvilinear_coordinate(line)
        step = 5 * np.mean(la_l[1:] - la_l[:-1])
        line = vmtk_resample_centerline(line, step)
        line = compute_splined_centerline(line, nknots=10, isline=True)
        area, sections = vmtk_compute_centerline_sections(capped_surface, line)

        # Get arrays
        a = get_point_data_array("CenterlineSectionArea", area)
        n = get_point_data_array("FrenetTangent", area, k=3)

        # Compute 'derivative' of the area
        dAdX = np.gradient(a.T[0], step)

        # Find the largest change in cross-sectional area
        tol = 5
        lim = -50
        ID = -1

        stop_id = np.nonzero(dAdX < lim)[0][ID] + tol
        normal = n[stop_id]
        center = area.GetPoint(stop_id)

        # Clip the model
        plane = vtk_plane(center, normal)
        surface, clipped = vtk_clip_polydata(surface, plane)

        # Find part to keep
        surface = vtk_clean_polydata(surface)
        clipped = vtk_clean_polydata(clipped)
        p_boundary = np.array(la_centerline_i.GetPoint(la_centerline_i.GetNumberOfPoints() - 1))

        surf_loc = get_vtk_point_locator(surface)
        clip_loc = get_vtk_point_locator(clipped)
        id_surf = surf_loc.FindClosestPoint(p_boundary)
        id_clip = clip_loc.FindClosestPoint(p_boundary)
        p_surface = np.array(surface.GetPoint(id_surf))
        p_clipped = np.array(clipped.GetPoint(id_clip))
        dist_surface = np.linalg.norm(p_surface - p_boundary)
        dist_clipped = np.linalg.norm(p_clipped - p_boundary)
        if dist_surface < dist_clipped:
            surface, clipped = clipped, surface

        surface = attach_clipped_regions_to_surface(surface, clipped, center)


if __name__ == '__main__':
    extract_LA_and_LAA()

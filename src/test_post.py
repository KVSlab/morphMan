from common import *

original_surface = read_polydata("example/surface.vtp")
original_cl = new_cl = read_polydata("example/surface_usr_centerline.vtp")
surface = read_polydata("example/surface_area_smooth_stenosis_2.0misr_0.5.vtp")

prepare_surface_output(surface, new_cl, original_surface, original_cl)

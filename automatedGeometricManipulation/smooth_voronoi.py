from common import *
from argparse import ArgumentParser
from os import path, listdir
from subprocess import STDOUT, check_output
from IPython import embed

import operator
import sys
import math
import numpy.linalg as la
import matplotlib.pyplot as plt

# Local import
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *
from moveandmanipulatetools import *

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('-d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('-case', type=str, default=None, help="Choose case")
    parser.add_argument('-s', '--smooth', type=bool, default=False,
                        help="If the original voronoi diagram (surface) should be" + \
                        "smoothed before it is manipulated", metavar="smooth")
    parser.add_argument("-a", "--alpha", type=float, default=None,
                        help="Compression factor in vertical direction, ranging from -1.0 to 1.0")
    parser.add_argument("-b", "--beta", type=float, default=None,
                        help="Compression factor in horizontal direction, ranging from -1.0 to 1.0")

    args = parser.parse_args()

    return args.smooth, args.dir_path, args.case, args.alpha, args.beta


def smooth_line_and_move_voronoi(case, dirpath, name, smoothingfactor, iterations, direction):
    it = iterations
    sf = smoothingfactor
    dr = direction
    # Input filenames
    model_path = path.join(dirpath, name, "model.vtp")

    # Output names
    model_smoothed_path = path.join(dirpath, name, "model_smoothed.vtp")
    new_centerlines_path = path.join("movedvoronoi", "%s_smooth_%s_%s_%s_centerlines.vtp" % (case, sf,dr,it))
    new_centerlines_path_tmp = path.join("movedvoronoi", "%s_smooth_%s_%s_%s_centerlines_tmp.vtp" % (case, sf,dr, it))
    model_new_surface = path.join("movedvoronoi", "%s_smooth_%s_%s_%s.vtp" % (case, sf,dr,it))
    model_new_surface_tmp = path.join("movedvoronoi", "%s_smooth_%s_%s_%s_tmp.vtp" % (case, sf,dr,it))

    # Centerlines
    centerline_complete_path = path.join(dirpath, name, "centerline_complete.vtp")
    centerline_clipped_path = path.join(dirpath, name,  "centerline_clipped.vtp")
    centerline_clipped_part_path = path.join(dirpath, name,  "centerline_clipped_part.vtp")
    centerline_clipped_part1_path = path.join(dirpath, name,  "centerline_clipped_end_part.vtp")
    centerline_new_path = path.join(dirpath, name, "centerline_interpolated.vtp")
    carotid_siphon_path = path.join(dirpath, name, "carotid_siphon.vtp")

    # Voronoi diagrams
    voronoi_path = path.join(dirpath, name, "model_voronoi.vtp")
    voronoi_smoothed_path = path.join(dirpath, name,  "model_voronoi_smoothed.vtp")
    voronoi_clipped_path = path.join(dirpath, name, "model_voronoi_clipped.vtp")
    voronoi_clipped_part_path = path.join(dirpath, name, "model_voronoi_clipped_part.vtp")

    # Read and check model
    if not path.exists(model_path):
        RuntimeError("The given directory: %s did not contain the file: model.vtp" % dirpath)

    # Clean surface
    surface = ReadPolyData(model_path)
    surface = surface_cleaner(surface)
    surface = triangulate_surface(surface)

    #Get a capped and uncapped version of the surface
    open_surface = surface
    capped_surface = capp_surface(surface)

    # Get inlet and outlets
    inlet, outlets = get_centers(open_surface, dirpath)
    
    # Compute all centerlines
    centerlines_complete = compute_centerlines(inlet, outlets,
                                               centerline_complete_path,
                                               capped_surface, resampling=0.1)

    print("Compute voronoi diagram")
    voronoi = makeVoronoiDiagram(surface, voronoi_path)
    if not path.exists(voronoi_smoothed_path) and smooth:
        voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerlines_complete, 0.25)
        WritePolyData(voronoi_smoothed, voronoi_smoothed_path)
        surface_smoothed = create_new_surface(voronoi_smoothed)
        WritePolyData(surface_smoothed, model_smoothed_path)

    voronoi = SmoothClippedVoronoiDiagram(voronoi, centerlines_complete, 0.25)

    centerlines_in_order = sort_centerlines(centerlines_complete)

    # Check if case includs the opthalmic artery
    #eye = False
    #eye, clip_ID, centerlines_in_order, eyeline = find_ophthalmic_artery(centerlines_in_order, clipping_points)

    #if eye == True:
    #    print("Clipping opthamlic artery")
    #    patch_eye = clip_eyeline(eyeline, clipping_points[0], clip_ID)


    # Get Carotid Siphon
    length = 0.01
    carotid_siphon = ReadPolyData(carotid_siphon_path)
    #carotid_siphon = CenterlineResampling(carotid_siphon, length)
    # NOTE: Add opthalmic and AcomArtery

    # Clipp Voronoi diagram
    masked_voronoi = MaskVoronoiDiagram(voronoi, carotid_siphon)
    clipped_voronoi = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi)
    
    # Smooth Carotid siphon 
    factor = smoothingfactor
    it = 300
    smooth_carotid_siphon = CenterlineSmoothing(carotid_siphon, factor, it)

    # Get rest of artery
    locator = get_locator(carotid_siphon)
    pStart = carotid_siphon.GetPoint(0)
    pEnd = carotid_siphon.GetPoint(carotid_siphon.GetNumberOfPoints() - 1)
    clipping_points = [pStart, pEnd]
    div_points = np.asarray(clipping_points)
    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)
    clip_points = points

    IDbif  = locator.FindClosestPoint(pEnd)

    patch_cl = CreateParentArteryPatches(centerlines_complete , clip_points)
    end_lines = []
    for i in range(patch_cl.GetNumberOfCells()):
        tmp_line  = ExtractSingleLine(patch_cl, i)
        if tmp_line.GetNumberOfPoints() > 50:
            end_lines.append(tmp_line)

    centerlines_end = merge_data(end_lines)
    masked_voronoi_end = MaskVoronoiDiagram(voronoi, centerlines_end)
    clipped_voronoi_end = ExtractMaskedVoronoiPoints(voronoi, masked_voronoi_end)

    # Update clipped Voronoi diagram
    moved_clipped_voronoi = move_voro(clipped_voronoi, carotid_siphon, smooth_carotid_siphon, direction)

    # Combine Voronoi diagram
    newVoronoi = merge_data([moved_clipped_voronoi,clipped_voronoi_end])

    # Create new surface
    new_surface = create_new_surface(newVoronoi)


    WritePolyData(new_surface, model_new_surface)
    return new_surface



def move_voro(clip_voro, old_cl, new_cl, direction):
    centerline_loc = get_locator(old_cl)

    newDataSet = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    verts = vtk.vtkCellArray()

    IDend = old_cl.GetNumberOfPoints() - 1
    IDstart_smooth = int(0.9*old_cl.GetNumberOfPoints())


    for p in range(clip_voro.GetNumberOfPoints()):
        cl_id = centerline_loc.FindClosestPoint(clip_voro.GetPoint(p))
        p0 = np.asarray(old_cl.GetPoint(cl_id))
        p1 = np.asarray(new_cl.GetPoint(cl_id))
        if direction == "smooth":
            dx = p1 - p0
        else:
            dx = -(p1 - p0)

        if cl_id > IDstart_smooth:
            dx = dx*(IDend - cl_id) / float(IDend - IDstart_smooth)

        dist = dx
        points.InsertNextPoint(np.asarray(clip_voro.GetPoint(p)) + dist)
        verts.InsertNextCell(1)
        verts.InsertCellPoint(p)

    newDataSet.SetPoints(points)
    newDataSet.SetVerts(verts)
    newDataSet.GetPointData().AddArray(clip_voro.GetPointData().GetArray(radiusArrayName))
    newVoronoi = newDataSet
    return newVoronoi


if  __name__ == "__main__":
    smooth, basedir, case, alpha, beta = read_command_line()
    folders = sorted([folder for folder in listdir(basedir) if folder[:2] in ["P0"]])
    name = "surface"
    smoothingfactor = 1.7
    iterations = 50
    directions = ["smooth"]#["smooth", "extend"]
    for folder in [folders[4]]:#, folders[8]]: 
        for i in range(1):
            print("==== Working on case %s ====" % folder)
            case = path.join(basedir,folder)
            direction = directions[i]
            smooth_line_and_move_voronoi(folder, case, name, smoothingfactor, iterations, direction)




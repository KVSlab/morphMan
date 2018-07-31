from common import *
from clipvoronoidiagram import SmoothClippedVoronoiDiagram
from os import listdir, path
from argparse import ArgumentParser
import sys


def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--dir_path', '-d', type=str, default=".",
                        help="Path to the folder with all the cases",
                        metavar="PATH")

    parser.add_argument('--case', '-c', type=str, default=None,
                        help="Name of spesific cases",
                        metavar="PATH")
    args = parser.parse_args()

    return args.dir_path, args.case

def start_case(casepath, folder):
    """
    Compute centerlines from inlet to two outlets,
    centerline between two outlets, and centerlines through
    complete geometry. 
    Compute Voronoi diagram of model, smooth Voronoi diagram, 
    and compute smoothed surface.
    
    Args:
        casepath (str): Location of cases.
        folder (str): Case name.
    """
    if folder.endswith(".tar.gz"):
        folder_ = path.join(casepath, folder[:-14])
    else:
        folder_ = path.join(casepath, folder)

    model_path = path.join(folder_, "surface", "model.vtp")
    model_smoothed_path = path.join(folder_, "surface", "model_smoothed.vtp")
    voronoi_path = path.join(folder_, "surface", "voronoi.vtp")
    voronoi_smoothed_path = path.join(folder_, "surface", "voronoi_smoothed.vtp")
    centerlines_path = path.join(folder_, "surface", "centerline_complete.vtp")
    centerlines_usr_path = path.join(folder_, "surface", "model_usr_centerline.vtp")
    centerlines_usr_bif_path = path.join(folder_, "surface", "model_usr_centerline_bif.vtp")
   
    # Untar
    if not path.isdir(folder_):
        check_output("bash untar.sh %s %s" % (casepath, folder), stderr=STDOUT, shell=True)
        check_output("cd -", stderr=STDOUT, shell=True)

    # Convert case info to plain text
    if not path.exists(path.join(folder_, "info.txt")):
        csv_to_txt(folder_)
        pass

    # Compute centerlines and voronoi diagram
    if not path.exists(centerlines_usr_path):
        print "Creating model_usr_centerline.vtp.."
        make_centerline(model_path, centerlines_usr_path, resampling=True, length=0.1, store_points=False)
        

    if not path.exists(centerlines_usr_bif_path):
        print "Creating model_usr_centerline_bif.vtp.."
        make_centerline(model_path, centerlines_usr_bif_path, resampling=True, length=0.1, store_points=False)

    if not path.exists(centerlines_path):
        print "Creating centerline_complete.vtp.."
        centerlines = make_centerline(model_path, centerlines_path, smooth=False,
                                    resampling=False, store_points=True)
    else:
        centerlines = read_polydata(centerlines_path)

    if not path.exists(voronoi_path):
        model = read_polydata(model_path)
        voronoi = make_voronoi_diagram(model, voronoi_path)
    else:
        voronoi = read_polydata(voronoi_path)

    if not path.exists(voronoi_smoothed_path):
        voronoi_smoothed = SmoothClippedVoronoiDiagram(voronoi, centerlines, 0.25)
        write_polydata(voronoi_smoothed, voronoi_smoothed_path)
    else:
        voronoi_smoothed = read_polydata(voronoi_smoothed_path)

    if not path.exists(model_smoothed_path):
        model_smoothed = create_new_surface(voronoi_smoothed)
        write_polydata(model_smoothed, model_smoothed_path)


if __name__ == "__main__":
    dir_path ,case = read_command_line()

    if case is None:
        home = "/home/henrik/master_uio/cases/" 

        cases = listdir(home)
        cases = [home + s for s in cases]
        for folder in cases:
            start_case(dir_path, folder)
    else:
        start_case(dir_path, case)

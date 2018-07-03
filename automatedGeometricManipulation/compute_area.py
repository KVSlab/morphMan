from common import *
from argparse import ArgumentParser
from os import path, listdir
from subprocess import STDOUT, check_output
from IPython import embed
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import resample, argrelextrema
import scipy.interpolate as interp

import operator
import sys
import math
import numpy.linalg as la
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams['axes.linewidth'] = 2.0

# Local import
from moveandmanipulatetools import *
from patchandinterpolatecenterlines import *
from clipvoronoidiagram import *
from paralleltransportvoronoidiagram import *

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


def get_clipping_points(dirpath):
    particles = path.join(dirpath, "optimal_bog.particles")
    all_points = np.loadtxt(particles)
    clipping_points = all_points[2:] 
    return clipping_points

def reconst_and_move_voro(dirpath, newmodelpath, newmodelclpath, smooth, name, method, k,K, alpha=0.0, beta=0.0):
    case = dirpath[-5:]
    # Input filenames
    model_path = path.join(dirpath, name, "model.vtp")

    # Find endID from landmarking
    siphon_path = path.join(dirpath, "surface", "carotid_siphon.vtp")

    # Output names
    model_smoothed_path = path.join(dirpath, name, "model_smoothed.vtp")
    model_new_surface = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s.vtp" % (alpha,beta))
    model_new_surface_tmp = path.join(dirpath, name, "manipulated_model_alpha_%s_beta_%s_tmp.vtp" % (alpha, beta))
    print dirpath
    new_centerlines_path = path.join("./meshmodels", "%s_centerlines_alpha_%s_beta_%s.vtp" % (case, alpha, beta))
    model_new_surface = path.join("./meshmodels", "%s_alpha_%s_beta_%s.vtp" % (case, alpha,beta))

    if method == "angle_sd":
        new_centerlines_path = path.join("./modelsangle", newmodelpath)
        model_new_surface = path.join("./modelsangle", newmodelclpath)
        new_centerlines_path = path.join("./meshmodels", newmodelpath)
        model_new_surface = path.join("./meshmodels", newmodelclpath)
    else:
        new_centerlines_path = path.join("./modelscurv", newmodelpath)
        model_new_surface = path.join("./modelscurv", newmodelclpath)
    # Centerlines
    centerline_complete_path = path.join(dirpath, name, "centerline_complete.vtp")
    centerline_clipped_path = path.join(dirpath, name,  "centerline_clipped.vtp")
    centerline_clipped_part_path = path.join(dirpath, name,  "centerline_clipped_part.vtp")
    centerline_clipped_part1_path = path.join(dirpath, name,  "centerline_clipped_end_part.vtp")
    centerline_new_path = path.join(dirpath, name, "centerline_interpolated.vtp")

    # Voronoi diagrams
    voronoi_path = path.join(dirpath, name, "model_voronoi.vtp")
    voronoi_smoothed_path = path.join(dirpath, name,  "model_voronoi_smoothed.vtp")
    voronoi_clipped_path = path.join(dirpath, name, "model_voronoi_clipped.vtp")
    voronoi_clipped_part_path = path.join(dirpath, name, "model_voronoi_clipped_part.vtp")
    
    # Extract Clipping points
    clipping_points = get_clipping_points(dirpath)

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

    voronoi = voronoi if not smooth else ReadPolyData(voronoi_smoothed_path)

    # Set clipping points
    div_points = np.asarray(clipping_points)
    points = vtk.vtkPoints()
    for point in div_points:
        points.InsertNextPoint(point)
    clip_points = points


    # Sort centerlines 
    lines = []
    n = centerlines_complete.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines_complete, i))

    longest = [lines[0]]
    longest_length = get_curvilinear_coordinate(longest[0])

    for i in range(1,n):
        tmplong = get_curvilinear_coordinate(lines[i])
        if len(tmplong) > len(longest_length):
            longest_length = tmplong
            longest.insert(0, lines[i])
        else: 
            longest.append(lines[i])

    centerlines_in_order = merge_data(longest)


    # Check if case includs the opthalmic artery
    eye = False
    eye, clip_ID, centerlines_in_order, eyeline = find_ophthalmic_artery(centerlines_in_order, clipping_points)


    # New model procedure 
    if method == "angle_sd":
        new_centerlines_path = path.join("./modelsangle", newmodelclpath)
        new_model_path = path.join("./modelsangle", newmodelpath)
        new_centerlines_path = path.join("./meshmodels", newmodelclpath)
        new_model_path = path.join("./meshmodels", newmodelpath)
    else:
        new_centerlines_path = path.join("./modelscurv", newmodelclpath)
        new_model_path = path.join("./modelscurv", newmodelpath)

    new_centerlines_complete = ReadPolyData(new_centerlines_path)
    new_surface = ReadPolyData(new_model_path)
    new_surface = surface_cleaner(new_surface)
    new_surface = triangulate_surface(new_surface)

    # Sort centerlines 
    lines = []
    n = new_centerlines_complete.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(new_centerlines_complete, i))

    longest = [lines[0]]
    longest_length = get_curvilinear_coordinate(longest[0])

    for i in range(1,n):
        tmplong = get_curvilinear_coordinate(lines[i])
        if len(tmplong) > len(longest_length):
            longest_length = tmplong
            longest.insert(0, lines[i])
        else: 
            longest.append(lines[i])

    new_centerlines_in_order = merge_data(longest)


    # Check if case includs the opthalmic artery
    eye = False
    eye, clip_ID, new_centerlines_in_order, eyeline = find_ophthalmic_artery(new_centerlines_in_order, clipping_points)
    #eye, clip_ID, new_centerlines_in_order, eyeline = find_ophthalmic_artery(centerlines_complete, clipping_points)
    new_centerlines_in_order = new_centerlines_complete


    print("Extract siphon") 
    centerline_to_change = ReadPolyData(siphon_path)
    endID = centerline_to_change.GetNumberOfPoints() - 1
    locator = get_locator(centerline_to_change)
    
    # Extract new line
    locator = get_locator(new_centerlines_complete)
    pN = centerline_to_change.GetPoint(endID)
    
    startID_new  = locator.FindClosestPoint(centerline_to_change.GetPoint(0))
    endID_new  = locator.FindClosestPoint(centerline_to_change.GetPoint(endID))
    if j < 7 or j in [12,13]:
        id = 4
    else:
        id = 0

    id =4
    line = new_centerlines_complete
    tmp_line = ExtractSingleLine(line,id)
    locator = get_locator(tmp_line)
    endID = locator.FindClosestPoint(pN)  

    #if (endID * 1.5) > tmp_line.GetNumberOfPoints():
    #    id += 1
    #    endID = get_locator(ExtractSingleLine(line,id)).FindClosestPoint(pN)

    new_centerline_to_change = ExtractSingleLine(new_centerlines_complete, id, endID=endID)
    misr = get_array(radiusArrayName, new_centerlines_complete)

    # Copy over MISR values
    if id != 0:
        start =0
        for i in range(id):
            start += ExtractSingleLine(new_centerlines_complete, i).GetNumberOfPoints()
    else:
        start = 0
    N = new_centerline_to_change.GetNumberOfPoints()
    size = ExtractSingleLine(new_centerlines_complete, id).GetNumberOfPoints()

    radiusArrayNumpy = misr[start:start+size][::-1][:endID+1]
    #radiusArrayNumpy = gaussian_filter(misr[start:start+size][::-1][:endID+1],2)

    radiusArray = get_vtk_array(radiusArrayName, 1, N)

    for i in range(N):
        radiusArray.SetTuple1(i, radiusArrayNumpy[i])

    new_centerline_to_change.GetPointData().AddArray(radiusArray)
    embed()
    # Compute Area
    print("Computing original area")
    centerline_splined = splineCenterline(centerline_to_change, 20)
    length=0.5
    centerline_splined_resamp = CenterlineResampling(centerline_splined, length, filename=None)
    centerline_area, section_old = compute_centerline_sections(surface, centerline_splined_resamp)
    print("Computing new area")
    new_centerline_splined = splineCenterline(new_centerline_to_change, 20)
    length = la.norm(np.array(centerline_splined.GetPoint(10)) - np.array(centerline_splined.GetPoint(11)))
    length=0.5
    new_centerline_splined_resamp = CenterlineResampling(new_centerline_splined, length, filename=None)
    new_centerline_area, section_new = compute_centerline_sections(new_surface, new_centerline_splined_resamp)

    if k == 0:
        if method[0] == "c":
            way = "cp"
        else:
            way = "ap"

    else:
        if method[0] == "c":
            way = "cm"
        else:
            way = "am"

    sec_old_path = "centerlinecomparison/section_%s_old.vtp" % (case)
    sec_new_path = "centerlinecomparison/section_%s_%s.vtp" % (case, way)

    WritePolyData(section_old, sec_old_path)
    WritePolyData(section_new, sec_new_path)
    return

    # Extract area

    area = get_array("CenterlineSectionArea", centerline_area)
    newarea = get_array("CenterlineSectionArea", new_centerline_area)
    
    length1, area1 = get_stats(centerline_area, dirpath, centerlines_complete)
    length2, area2 = get_stats(new_centerline_area, dirpath, new_centerlines_complete)

    locator1 = get_locator(centerline_area)
    locator2 = get_locator(new_centerline_area)

    p1 = clipping_points[0]
    p2 = clipping_points[1]

    ID1 = locator1.FindClosestPoint(p1)
    ID2 = locator1.FindClosestPoint(p2)
    ID1_ = locator2.FindClosestPoint(p1)
    ID2_ = locator2.FindClosestPoint(p2)

    siphon_area = area1[ID1:ID2]
    newsiphon_area = area2[ID1_:ID2_]

    area_interp = interp.interp1d(np.arange(siphon_area.size),siphon_area.T)
    siphon_area2 = area_interp(np.linspace(0,siphon_area.size-1,newsiphon_area.size))

    # Compare with MISR area
    misr_old = get_array(radiusArrayName, centerline_splined)
    misr_new = get_array(radiusArrayName, new_centerline_splined_resamp)
    misr_old_area = np.pi * misr_old ** 2
    misr_new_area = np.pi * misr_new ** 2

    siphon_misr_area = misr_old_area[ID1:ID2]
    newsiphon_misr_area = misr_new_area[ID1_:ID2_]

    misr_area_interp = interp.interp1d(np.arange(siphon_misr_area.size),siphon_misr_area.T)
    siphon_misr_area2 = misr_area_interp(np.linspace(0,siphon_misr_area.size-1,newsiphon_misr_area.size))

    plt.figure(figsize=(10,5))

    plt.plot(siphon_misr_area2.T, 'r-', label="Original area (MISR)",linewidth=3)
    plt.plot(newsiphon_misr_area, '-', label="New area (MISR)", color="gray", linewidth=3)

    plt.plot(siphon_area2.T, 'r--', label="Original area",linewidth=3)
    plt.plot(newsiphon_area, '--', label="New area", color="grey", linewidth=3)
    #plt.margins(0)
    plt.ylabel(r"Area [mm$^2$]", rotation=90 ,fontsize=25, labelpad =20)
    plt.xlabel(r"Abscissa, Anterior bend",fontsize=25, labelpad=20)
    plt.yticks(fontsize=10)#, rotation = 270)
    plt.xticks(fontsize=10)#, rotation = 270)
    plt.legend()
    #plt.show()

    if k == 0:
        if method[0] == "c":
            way = "cp"
            way = "ap"
        else:
            way = "ap"

    else:
        if method[0] == "c":
            way = "cm"
            way = "ap"
        else:
            way = "am"

    sec_old_path = "centerlinecomparison/section_%s_old.vtp" % (case)
    sec_new_path = "centerlinecomparison/section_%s_%s.vtp" % (case, way)

    WritePolyData(section_old, sec_old_path)
    WritePolyData(section_new, sec_new_path)

    plt.savefig("centerlinecomparison/misr/%s_%s_%s.png"  %(case,method,way),dpi=300, format='png', bbox_inches="tight") 
    plt.close("all")


def CenterlineResampling(line, length, filename=None):

    if filename is None:
        filename = "tmp_re.vtp"
        WritePolyData(line, filename)

    command = ('vmtk vmtkcenterlineresampling -ifile %s -ofile %s -length %s') % (filename, filename, length)

    a = check_output(command, stderr=STDOUT, shell=True)

    status, text = success(a)

    if not status:
        print "Something went wrong when resampling the centerline"
        print text
        sys.exit(1)

    line = ReadPolyData(filename)
    
    return line 


def cutLine(line, endID):
    # Create edges between new_centerline points
    pts = vtk.vtkPoints()
    for i in range(endID):
        pts.InsertNextPoint(line.GetPoint(i))

    lines = vtk.vtkCellArray()
    for i in range(endID-2):
        newline = vtk.vtkLine()
        newline.GetPointIds().SetId(0, i)
        newline.GetPointIds().SetId(1, i + 1)
        lines.InsertNextCell(newline)

    line = vtk.vtkPolyData()
    line.SetPoints(pts)
    line.SetLines(lines)
    return line


def find_ophthalmic_artery(centerlines, clip_pts):

    # Extract lines:
    lines = []
    n = centerlines.GetNumberOfCells()
    for i in range(n):
        lines.append(ExtractSingleLine(centerlines, i))
    #    WritePolyData(lines[i], "line%i.vtp"%i)

    longest = lines[0]
    tol = 0.40
    
    # Find start and stop IDs along clipped curve
    locator_longest = get_locator(longest)
    p1      = clip_pts[0]
    p2      = clip_pts[1]
    ID1     = locator_longest.FindClosestPoint(p1)
    ID2     = locator_longest.FindClosestPoint(p2)
    
    eye = False
    index = 1
    for line in lines[1:]:
        locator_checkline = get_locator(line)
        len_check = len(get_curvilinear_coordinate(line)) 
        if len_check < ID2:
            IDStop = len_check - 1
        else:
            IDStop = ID2

        for i in np.arange(ID1, IDStop):
            p_eye   = np.asarray(line.GetPoint(i))
            p_cl    = np.asarray(longest.GetPoint(i))
            dist = la.norm(p_eye - p_cl)
            if  dist > tol:
                clip_ID = i
                centerlines = merge_data(lines)
                eye = True
                eyeline = lines[index]
                del lines[index]
                centerlines = merge_data(lines)
                break
        if eye:
            break
        index += 1
    
    if eye:
        return eye, clip_ID, centerlines, eyeline
    else:
        return eye, None, centerlines, None

   

def get_lineToChange(centerline, tol, start=0):
    line2 = ExtractSingleLine(centerline, start)
    numberOfPoints2 = line2.GetNumberOfPoints()

    n = start + 2
    pointIDs = []
    for j in range(start+1, n):
        line1 = ExtractSingleLine(centerline, j)
        numberOfPoints1 = line1.GetNumberOfPoints()

        N = min(numberOfPoints1, numberOfPoints2)
        for i in range(N):
            point1 = line1.GetPoints().GetPoint(i)
            point2 = line2.GetPoints().GetPoint(i)
            if math.sqrt(distance(point1, point2)) > tol:
                pointID = i
                break

        pointIDs.append(pointID) #move_past_sphere(line2, pointID))

    pointID = min(pointIDs)
    lineToChange = ExtractSingleLine(centerline, start, endID=pointID)
    lineToChange = splineCenterline(lineToChange, knots=25)
    
    #curvature = get_array("Curvature", line)
    #curvature = get_array("Torsion", line)
    #length = get_curvilinear_coordinate(line)
    
    return lineToChange

def splineCenterline(line, knots):
    # Allocate data structure to store centerline points
    data = np.zeros((line.GetNumberOfPoints(), 3))
    MISR = get_array(radiusArrayName, line)

    # Collect data from centerline
    curv_coor = get_curvilinear_coordinate(line)
    for i in range(data.shape[0]):
        data[i,:3] = line.GetPoints().GetPoint(i)

    t = np.linspace(curv_coor[0], curv_coor[-1], knots+2)[1:-1]
    fx = splrep(curv_coor, data[:,0], k=4, t=t)
    fy = splrep(curv_coor, data[:,1], k=4, t=t)
    fz = splrep(curv_coor, data[:,2], k=4, t=t)

    fx_ = splev(curv_coor, fx)
    fy_ = splev(curv_coor, fy)
    fz_ = splev(curv_coor, fz)

    # Store data for converting to vtkPolyLine
    data = np.zeros((len(curv_coor), 4))
    data[:,0] = fx_
    data[:,1] = fy_
    data[:,2] = fz_
    data[:,3] = MISR[:,0]

    header = ["X", "Y", "Z", radiusArrayName]
    line = data_to_vtkPolyData(data, header)

    return line


def get_stats(centerline_area, folder, centerline):
    area = get_array("CenterlineSectionArea", centerline_area)
    MISR_ = get_array(radiusArrayName, centerline)**2*math.pi
    MISR = np.array([MISR_[i] for i in range(MISR_.shape[0] - 1, -1, -1)])[:area.shape[0]]
    length = get_curvilinear_coordinate(centerline_area)

    for i in range(2):
        area = gaussian_filter(area, 5)
        MISR = gaussian_filter(MISR, 5)
    
    # Area shape - "circleness"
    circleness = MISR / area
    max_circleness = circleness.max()
    min_circleness = circleness.min()
    mean_circleness = circleness.mean()
        
    # Local extremas, ignore min or max on boundary
    local_min_MISR_ID = argrelextrema(MISR, np.less)[0]
    local_max_MISR_ID = argrelextrema(MISR, np.greater)[0]
    local_min_area_ID = argrelextrema(area, np.less)[0]
    local_max_area_ID = argrelextrema(area, np.greater)[0]

    local_min_MISR = MISR[local_min_MISR_ID]
    local_max_MISR = MISR[local_max_MISR_ID]
    local_min_area = area[local_min_area_ID]
    local_max_area = area[local_max_area_ID]

    global_min_MISR = local_min_MISR.min()
    global_max_MISR = local_max_MISR.max()
    global_min_area = local_min_area.min()
    global_max_area = local_max_area.max()

    mean_area = area.mean()
    mean_MISR = MISR.mean()

    number_of_max = local_max_area.shape[0]

    # Min max derived parameters
    length_min_max = abs(length[(area == global_min_area).nonzero()[0]] - \
                         length[(area == global_max_area).nonzero()[0]])[0]

    max_mean_ratio_area = global_max_area / mean_area
    min_mean_ratio_area = global_min_area / mean_area
    max_mean_ratio_MSIR = global_max_MISR / mean_MISR
    min_mean_ratio_MISR = global_min_MISR / mean_MISR

    # Global and local disent
    global_min_max_disent = abs(math.sqrt(global_max_area)/math.pi -
                                math.sqrt(global_min_area) / math.pi) / \
                            length_min_max
    local_max_stepest_disent = 0

    if length[(area == local_min_area[0]).nonzero()[0]] > \
       length[(area == local_max_area[0]).nonzero()[0]]:
        start = 1
    else:
        start = 0

    N = min(number_of_max, local_min_area.shape[0] - start)
    for i in range(N):
        min_ = local_min_area[start + i]
        max_ = local_max_area[i]

        h = math.sqrt(max_)/math.pi - math.sqrt(min_) / math.pi
        l = abs(length[(area == max_).nonzero()[0]] - length[(area == min_).nonzero()[0]])
        if h/l > local_max_stepest_disent:
            local_max_stepest_disent = h / l
    
    
    # Max point disent (max |derivative|)
    knots = 40
    t = np.linspace(length[0], length[-1], knots+2)[1:-1]
    spline_area = splrep(length, area, k=4, t=t)
    spline_area_ = splev(length, spline_area)
    darea_dx = splev(length, spline_area, der=1)
    max_derivative = abs(darea_dx).max()

    #print "max_min_ratio_area:", global_max_area / global_min_area
    stats = {"max_derivative": max_derivative,
             "local_max_stepest_disent": local_max_stepest_disent,
             "max_mean_ratio_area": max_mean_ratio_area,
             "min_mean_ratio_area" : min_mean_ratio_area,
             "mean_area": mean_area,
             "max_min_ratio_area": global_max_area / global_min_area,
             "length_min_max": length_min_max,
             "global_min_area": global_min_area,
             "global_max_area": global_max_area,
             "max_circleness": max_circleness,
             "min_circleness": min_circleness,
             "mean_circleness": mean_circleness,
             "number_of_max": number_of_max}

    writeParameters(stats, folder)
    
    return length, area


if  __name__ == "__main__":
    smooth, basedir, case, alpha, beta = read_command_line()
    folders = sorted([folder for folder in listdir(basedir) if folder[:2] in ["P0"]])
    name = "surface"
    growDir = "both"
    rad = 0.2

    a_sd = "plots_angle/optimal_alphabeta.txt"
    model_a_sd = "modelsangle/"
    model_a_sd = "meshmodels/"

    c_sd = "plots_curv/optimal_alphabeta.txt"
    model_c_sd = "modelscurv/"
    a1 = 0.89449081803
    b1=-0.188013355593
    a2=-0.2
    b2=0.638430717863
    a_ = [a1,a2]
    b_ = [b1,b2]
    cases =sorted(listdir("cases/"))
    methods = ["angle_sd", "curv_sd"]#, "angle_05sd", "curv_sd", "curv_2sd"]
    models = [model_a_sd, model_c_sd]#, model_a_sd_half, model_c_sd, model_c_sd_double]
    params = [a_sd, c_sd]#, a_sd_half, c_sd, c_sd_double]

    
    for i,method in enumerate(methods[0:1]):
        files = sorted([f for f in listdir(models[i]) if "centerlines" not in f])
        k = 0
        j = 0
        for jj,f in enumerate(files[j:]):
            
            folder = f[:5]
            info = f.split("_")
            alpha = float(info[2])
            beta = float(info[-1][:-4])
            newmodelpath = f
            newmodelclpath = f[:6] + "centerlines_" + f[6:]
            print alpha,beta
            
            print("==== Working on case %s ====" % folder)
            case = path.join(basedir,folder)
            reconst_and_move_voro(case ,newmodelpath, newmodelclpath, smooth, name, method, k, j, alpha, beta)
            if k == 0: k = 1
            else: k = 0
            j+=1
            if jj > 1:
                break




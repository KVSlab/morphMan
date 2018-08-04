import xml.etree.ElementTree as ET
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

from matplotlib import rc, rcParams
from os import remove
from common import *


paraview.simple._DisableFirstRenderCameraReset()
rc('text', usetex=True)
rc('font', family='serif')
rcParams['axes.linewidth'] = 2.0
rcParams['figure.figsize'] = 9.2, 12

def set_camera_position(renderView1, filename):
    """
    Adjust camera position to given
    position saved in input file .pvcc
    from Paraview
    """
    tree = ET.parse(filename)
    root = tree.getroot()

    # FIXME: Not sure if index is sorted, and thus the need for filling with 0
    position = [0, 0, 0]
    view_up = [0, 0, 0]
    focal = [0, 0, 0]
    rotation = [0, 0, 0]

    for child in root[0]:
        if child.attrib['name'] == 'CameraFocalPoint':
            for subChild in child:
                focal[int(subChild.attrib["index"])] = float(subChild.attrib['value'])
        if child.attrib['name'] == 'CenterOfRotation':
            for subChild in child:
                rotation[int(subChild.attrib["index"])] = float(subChild.attrib['value'])
        if child.attrib['name'] == 'CameraPosition':
            for subChild in child:
                position[int(subChild.attrib["index"])] = float(subChild.attrib['value'])
        if child.attrib['name'] == 'CameraViewUp':
            for subChild in child:
                view_up[int(subChild.attrib['index'])] = float(subChild.attrib['value'])
        if child.attrib["name"] == "CameraParallelScale":
            for subChild in child:
                parallel_scale = float(subChild.attrib['value'])
        if child.attrib["name"] == "RotationFactor":
            for subChild in child:
                rotation_factor = float(subChild.attrib['value'])
        if child.attrib["name"] == "CameraViewAngle":
            for subChild in child:
                view_angle = float(subChild.attrib['value'])
        if child.attrib["name"] == "CameraParallelProjection":
            for subChild in child:
                if "value" in subChild.attrib.keys():
                    parallel_projection = bool(subChild.attrib['value'])

    renderView1.CenterOfRotation = rotation
    renderView1.CameraPosition = position
    renderView1.CameraFocalPoint = focal
    renderView1.CameraViewUp = view_up
    renderView1.CameraParallelScale = parallel_scale
    renderView1.CameraParallelProjection = 0 #parallel_projection
    renderView1.CameraViewAngle = view_angle
    renderView1.RotationFactor = rotation_factor
    renderView1.CameraViewUp = view_up

    return renderView1

def combine_landmark(method):
    """
    Combine landmarking results into an
    overall comparison image.

    Args:
        method (str): Landmarking method
    """
    savepath = "/home/henrik/article/results/"
    savedir = "../results/landmark_%s.png" % method

    cases = sorted(listdir(savepath))
    images = [path.join(savepath, img) for img in sorted(listdir(savepath))]
    n_images = len(images)

    fig = plt.figure()
    cols = int(np.ceil(len(cases)/2.))
    for n, image in enumerate(images):
        a = fig.add_subplot(cols, np.ceil(n_images/float(cols)), n + 1)
        img = mpimg.imread(image)
        plt.imshow(img)
        plt.xticks([])
        plt.yticks([])
        plt.ylabel(r"%s" % cases[n][:-4],fontsize=39)
    plt.tight_layout()
    plt.savefig(savedir, format="png", bbox_inches="tight", dpi=300)


def viz_landmark(case, points_name):
    """
    Visualize and save screenshot of models
    including centerline and landmarking interfaces,
    visualized in Paraview.

    Args:
        case (str): Case name.
        points_name (str): Name of file with landmarking points.

    """

    home = "/home/henrik/article"
    pointpath = "%s/cases/%s/%s" % (home, case, points_name)
    siphonpath = "%s/cases/%s/surface/carotid_siphon.vtp" % (home, case)
    modelpath = "%s/cases/%s/surface/model.vtp" % (home, case)
    savepath = "../results/%s.png" % ( case)

    # Adjust points to lie on centerline
    try:
        tmp_particles = "tmp.particles"
        siphon = read_polydata(siphonpath)
        output = open(pointpath, "r")
        lines = output.readlines()
        P_old = []
        P_new = []
        for line in lines:
            p = line.split()
            P = np.array([float(p[0]),float(p[1]),float(p[2])])
            P_old.append(P)

        locator = get_locator(siphon)
        for p in P_old:
            P = siphon.GetPoint(locator.FindClosestPoint(p))
            P_new.append(P)

        output = open(tmp_particles, "w")
        for p in P_new:
            output.write("%s %s %s\n" % (p[0], p[1], p[2]))

        output.close()
    except IOError:
        pass

    # create a new 'XML PolyData Reader'
    carotid_siphonvtp = XMLPolyDataReader(FileName=[siphonpath])
    modelvtp = XMLPolyDataReader(FileName=[modelpath])

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    renderView1.ViewSize = [1442, 860]

    # show data in view
    modelvtpDisplay = Show(modelvtp, renderView1)
    # trace defaults for the display properties.
    modelvtpDisplay.ColorArrayName = [None, '']
    modelvtpDisplay.GlyphType = 'Arrow'

    # reset view to fit data
    renderView1.ResetCamera()
    renderView1.OrientationAxesVisibility = 0

    # show data in view
    carotid_siphonvtpDisplay = Show(carotid_siphonvtp, renderView1)
    # trace defaults for the display properties.
    carotid_siphonvtpDisplay.ColorArrayName = [None, '']
    carotid_siphonvtpDisplay.GlyphType = 'Arrow'

    # Properties modified on modelvtpDisplay
    modelvtpDisplay.Opacity = 0.4

    # Properties modified on modelvtpDisplay
    modelvtpDisplay.LineWidth = 10.0

    # set active source
    SetActiveSource(carotid_siphonvtp)

    # Properties modified on carotid_siphonvtpDisplay
    carotid_siphonvtpDisplay.LineWidth = 10.0

    # change solid color
    carotid_siphonvtpDisplay.DiffuseColor = [1.0, 0.0, 0.01568627450980392]

    # create a new 'Particles Reader'
    optimal_bogparticles = ParticlesReader(FileName=[tmp_particles])

    # get color transfer function/color map for 'Scalar'
    scalarLUT = GetColorTransferFunction('Scalar')

    # show data in view
    optimal_bogparticlesDisplay = Show(optimal_bogparticles, renderView1)
    # trace defaults for the display properties.
    optimal_bogparticlesDisplay.ColorArrayName = ['POINTS', 'Scalar']
    optimal_bogparticlesDisplay.LookupTable = scalarLUT
    optimal_bogparticlesDisplay.GlyphType = 'Arrow'

    # show color bar/color legend
    optimal_bogparticlesDisplay.SetScalarBarVisibility(renderView1, False)

    # get opacity transfer function/opacity map for 'Scalar'
    scalarPWF = GetOpacityTransferFunction('Scalar')

    # change representation type
    optimal_bogparticlesDisplay.SetRepresentationType('3D Glyphs')

    # Properties modified on optimal_bogparticlesDisplay
    optimal_bogparticlesDisplay.GlyphType = 'Sphere'

    # turn off scalar coloring
    ColorBy(optimal_bogparticlesDisplay, None)

    # change solid color
    optimal_bogparticlesDisplay.DiffuseColor = [0.0, 0.0, 0.0]

    # Properties modified on optimal_bogparticlesDisplay.GlyphType
    optimal_bogparticlesDisplay.GlyphType.Radius = 0.5

    # Properties modified on renderView1
    renderView1.Background = [1.0, 1.0, 1.0]

    #### saving camera placements for all active views

    # current camera placement for renderView1

    # reset view to fit data
    camerapath = "../camera/%s.pvcc" % case[2:]
    set_camera_position(renderView1,camerapath)
    paraview.simple._DisableFirstRenderCameraReset()
    #### uncomment the following to render all views
    #RenderAllViews()
    SaveScreenshot(savepath, magnification=1, quality=500)
    # alternatively, if you want to write images, you can use SaveScreenshot(...).


    Delete(carotid_siphonvtp)
    Delete(modelvtp )
    Delete(optimal_bogparticles )
    Delete(carotid_siphonvtpDisplay)
    Delete(modelvtpDisplay)
    Delete(optimal_bogparticlesDisplay)
    del carotid_siphonvtp
    del modelvtp
    del optimal_bogparticles
    del carotid_siphonvtpDisplay
    del modelvtpDisplay
    del optimal_bogparticlesDisplay

    try:
        remove(tmp_particles)
    except OSError:
        pass

def main():
    cases = sorted(listdir("../cases"))
    points = "carotid_siphon_points.particles"

    method = "bogunovic"
    for case in cases:
        viz_landmark(case, points)

    combine_landmark(method)

if __name__ == "__main__":
    main()

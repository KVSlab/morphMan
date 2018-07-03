from paraview.simple import *
from os import path, listdir,  rename, makedirs
from IPython import embed
import xml.etree.ElementTree as ET
import re

import numpy as np
import sys
import math
import os
import numpy as np
import matplotlib.pyplot as plt

from os import path, sep, listdir
from hashlib import sha1
from argparse import ArgumentParser
from matplotlib import rc, rcParams
import matplotlib.image as mpimg
rc('text', usetex=True)
rc('font', family='serif')
rcParams['axes.linewidth'] = 1.0
rcParams['figure.figsize'] = 8, 10



def combine_plot(case, method):
    pass

def set_camera_position(renderView1, filename):
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


def plot(oldmodel,newmodel, k, method):
    methodpath = "%smodels" % method
    initial_path = path.join("initialmodels", oldmodel)
    moved_path = path.join(methodpath, newmodel)
    print initial_path
    print moved_path
    case = newmodel.split("_")[0]

    sign = newmodel.split("_")[-1].split(".")[0][0]

    if sign == "-":
        if method[0] == "c":
            way = "plus"
        else:
            way = "minus"


    else:
        if method[0] == "c":
            way = "minus"
        else:
            way = "plus"

    savepath = "curvmodelscompare/%s_%s_%s.png" % (case, method, way)
    savepath = "anglemodelscompare/%s_%s_%s.png" % (case, method, way)
    # create a new 'XML PolyData Reader'
    moved = XMLPolyDataReader(FileName=[moved_path])
    old = XMLPolyDataReader(FileName=[initial_path])
    moved.PointArrayStatus = ['Scalars_', 'Normals_']


    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')
    # uncomment following to set a specific view size
    renderView1.ViewSize = [1214, 860]
    renderView1.OrientationAxesVisibility = 0

    camerapath = "../results/camera/compare/%s.pvcc" % case[2:]
    set_camera_position(renderView1,camerapath)

    # get color transfer function/color map for 'Scalars'
    scalarsLUT = GetColorTransferFunction('Scalars')

    # show data in view
    movedDisp = Show(moved, renderView1)
    # trace defaults for the display properties.
    movedDisp.ColorArrayName = ['POINTS', 'Scalars_']
    movedDisp.LookupTable = scalarsLUT
    movedDisp.GlyphType = 'Arrow'


    # show color bar/color legend
    movedDisp.SetScalarBarVisibility(renderView1, False)

    # get opacity transfer function/opacity map for 'Scalars'
    scalarsPWF = GetOpacityTransferFunction('Scalars')

    # create a new 'XML PolyData Reader'

    #old.CellArrayStatus = ['Array 0x160887b50']

    # show data in view
    oldDisp = Show(old, renderView1)
    # trace defaults for the display properties.
    oldDisp.ColorArrayName = [None, '']
    oldDisp.GlyphType = 'Arrow'

    # set active source
    SetActiveSource(moved)

    # turn off scalar coloring
    ColorBy(movedDisp, None)

    # Properties modified on renderView1
    renderView1.Background = [1.0, 1.0, 1.0]

    # change solid color
    movedDisp.DiffuseColor = [1.0, 0.0, 0.01568627450980392]


    # set active source
    SetActiveSource(old)

    # Properties modified on model_P0086vtp_1Display
    oldDisp.Opacity = 0.5

    # set active source
    SetActiveSource(moved)

    # Properties modified on p0086_alpha_01485_beta_00348vtpDisplay
    movedDisp.Opacity = 0.5

    #### saving camera placements for all active views

    # current camera placement for renderView1
    # Hide orientation axes

    #### uncomment the following to render all views
    RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    SaveScreenshot(savepath)


    Delete(old)
    Delete(oldDisp)
    Delete(moved)
    Delete(movedDisp)
    del old
    del oldDisp
    del moved
    del movedDisp


def plot_all_features():
    curv = "curvmodelscompare"
    angle = "anglemodelscompare"
    folders_curv = sorted(listdir(curv))
    folders_angle = sorted(listdir(angle))

    #images = [path.join(savepath, img) for img in sorted(listdir(savepath)) if "Wpng" in img]
    images_curvp = [path.join(curv, s) for s in folders_curv if "plus" in s]
    images_curvm = [path.join(curv, s) for s in folders_curv if "minus" in s]
    images_anglep = [path.join(angle, s) for s in folders_angle if "plus" in s]
    images_anglem = [path.join(angle, s) for s in folders_angle if "minus" in s]

    meth = [curv, angle]
    cases = sorted(listdir("cases"))
    ways = ["plus", "plus", "minus", "minus"]
    for i,images in enumerate([images_curvp, images_anglep, images_curvm, images_anglem]):
        savedir = "%s_%s.png" % (meth[i % 2], ways[i])

        n_images = len(images)
        fig = plt.figure()
        cols = 5
        for n, image in enumerate(images):
            ax = fig.add_subplot(cols, np.ceil(n_images/float(cols)), n + 1)
            img = mpimg.imread(image)
            plt.imshow(img)
            plt.xticks([])
            plt.yticks([])
            plt.ylabel(r"%s" % cases[n],fontsize=28, labelpad =5)
            #if n < 3:
            #    plt.xlabel(r"%s" % ways[n],fontsize=23, labelpad =10)
            #    ax.xaxis.set_label_position('top') 
            #    ax.xaxis.tick_top()

            plt.subplots_adjust(wspace = -0.21, hspace = 0.10)

        #plt.suptitle(r"%s" % names[i], fontsize=31, y = 0.95, x=0.513)
        #plt.tight_layout()
        plt.savefig(savedir, format="png", bbox_inches="tight", dpi=300)
        #plt.show()




def main():
    RenderAllViews()
    curvcases =sorted(listdir("./curvmodels/"))
    anglecases =sorted(listdir("./anglemodels/"))
    cases =sorted(listdir("./initialmodels/"))
    #for i in range(2): #Method iteration
    i = 0
    method = "angle"
    #method = "curv"
    for j, case in enumerate(cases): # Case iteration
        for k in range(2):  # plus / minus iteration
            newmodel = anglecases[i]
            #newmodel = curvcases[i]
            oldmodel = cases[j]
            i += 1
           # plot(oldmodel,newmodel, k, method)

    plot_all_features()

main()

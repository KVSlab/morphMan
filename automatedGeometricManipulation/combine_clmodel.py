import Image

from os import path, listdir,  rename, makedirs
from IPython import embed

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import rc, rcParams
from scipy import ndimage
from matplotlib import gridspec

# plot it

rc('text', usetex=True)
rc('font', family='serif')
rcParams['axes.linewidth'] = 0.0
rcParams['figure.figsize'] = 15, 5

def combine_plot(case, method,m):
    models = "/home/henrik/master_uio/manipulate/powerpoint/modelscompare"
    cldir = "/home/henrik/master_uio/manipulate/centerlinecomparison"
    
   # cases = sorted(listdir(models))
    imagesmodel = [path.join(models, img) for img in sorted(listdir(models)) if "angle" in img][m:m+2]
    imagescl = [path.join(cldir, img) for img in sorted(listdir(cldir)) if "angle" in img][m:m+2]
    images = [imagescl[0], imagesmodel[0], imagesmodel[1], imagescl[1]]
    if m > 10:
        embed()
    n_images = 4
    cols = 1
    #print images
    rect_box = [0,40 , 0, 40]
    fig = plt.figure()
    for n, image in enumerate(images):

        a = fig.add_subplot(cols, np.ceil(n_images/float(cols)), n + 1)
        img = mpimg.imread(image)
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        if n in[1,2]:
            plt.imshow(img, aspect="auto")
        else:
            plt.imshow(img,aspect=1.0)

        plt.subplots_adjust(wspace=0.0)

        #plt.ylabel("%s" % cases[n][:-4],fontsize=39)
    fig.set_size_inches(np.array(fig.get_size_inches()) * n_images)
    #plt.show()
    savedir = "clandmodels/%s_%s.png" % (case,method)
    plt.savefig(savedir, format="png", bbox_inches="tight", dpi=100)



def main():

    cases =sorted(listdir("cases/"))
    method = ["curv_sd","angle_sd"]#, "angle_05sd", "curv_sd", "curv_2sd"]
    m = 16
    for case in cases[8:]:
        combine_plot(case, method[1],m)
        m+=2


main()

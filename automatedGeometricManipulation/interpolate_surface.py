# -*- coding: utf-8 -*-
from pylab import rc, rcParams
from os import listdir, path
from scipy.interpolate import griddata
from scipy import interpolate
from IPython import embed

import numpy as np
import numpy.linalg as la

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

rcParams['figure.figsize'] = 15, 10

# Find intersection
def intersect(method, f,alpha_long,  beta_long,tol=None):
    zeros = []
    for i in alpha_long:
        for j in beta_long:
            diff = abs(f(j,i) - method(0,0, tol))
            if "c" in method.__name__:
                if diff < 0.001:
                    zeros.append([i,j])
            else:
                if diff < 0.05:
                    zeros.append([i,j])

    
    return zeros

def optimal_curvature(i, files, plot):
    files = path.join("param_curv", files)
    method = "curv"
    alphabeta_bound = np.loadtxt("alphabeta_bound.txt")

    #with open('new_allcases_initialvalues.txt') as file:
    #    initval = [[float(digit) for digit in line.split()] for line in file]

    amin, amax, bmin, bmax = alphabeta_bound[i][0],alphabeta_bound[i][1],alphabeta_bound[i][2],alphabeta_bound[i][3] 

    #curv0 = initval[i][0]
    sd_curv = 0.045

    # Defined SD planes for curvature
    # Tolerance added for adjusting SD 
    # if there are no intersections found
    def cpsd(x, y, tol=None):
        if tol is not None:
            return curv0 + sd_curv - tol 
        else:
            return curv0 + sd_curv 

    def cmsd(x, y, tol=None):
        if tol is not None:
            return curv0 - sd_curv + tol 
        else:
            return curv0 - sd_curv 
        
    def c(x, y, tol=None):
        return curv0

    N = 50
    edge = [amin,amax,bmin,bmax]
    alpha = np.linspace(amin,amax,N)
    beta = np.linspace(bmin,bmax,N) 
    alpha_long = np.linspace(amin,amax,300)
    beta_long = np.linspace(bmin,bmax,300) 
    
    with open(files) as file:
        curv = [[float(digit) for digit in line.split()] for line in file]

    yy,xx = np.meshgrid(beta,alpha)

    points = np.zeros((N,2))
    for i in range(len(xx)):
        points[i] = [alpha[i],beta[i]]

    f = interpolate.interp2d(beta,alpha, curv, kind='cubic')
    curv0 = f(0,0)
    zz = f(beta, alpha)

    case = files[-9:-4]
    print "Working on %s" % case
    Z = []
    for plane in [cpsd, cmsd, c]: 
        print "Method: %s" % (plane.__name__)
        # Find intersecting points
        # Reduces SD if no points are found
        zeros = intersect(plane, f, alpha_long,beta_long)
        if len(zeros) > 10:
            dx = int(len(zeros) / 10.)
        elif len(zeros) > 0:
            dx = 1
        else:
            empty = True
            tol = 0.005
            maxiter = 50
            it = 0
            print "Found no points..Adjusting SD"
            while empty == True and it < maxiter:
                print "Iterations: %i" % (it + 1)
                zeros = intersect(plane,f, alpha_long, beta_long, tol)
                if len(zeros) > 0:
                    empty = False
                it += 1
                tol += 0.001

            if len(zeros) > 10:
                dx = int(len(zeros) / 10.)
            elif len(zeros) > 0:
                dx = 1

        if len(zeros) > 0:

            points = []
            for p in zeros:
                if plane.__name__ == "cpsd":
                    if p[1] < 0:
                        points.append(p)
                elif plane.__name__ == "cmsd":
                    if p[1] > 0:
                        points.append(p)

            if plane != c:
                P = points[0]
                tol = 0.15
                dist = 10
                for p in points[1:]:
                    dist_tmp = la.norm(np.array(p))
                    if dist_tmp < dist and dist_tmp > tol:
                        dist = dist_tmp
                        P = p
                Z.append(P)
                writeOptimalP(case, P, plane.__name__)

            #plotAlphaBeta(zeros, case, plane.__name__, edge)
            
    plus = f(Z[0][1], Z[0][0]) [0]
    minus = f(Z[1][1], Z[1][0]) [0]
    init = curv0[0]
    #return plus,minus,init

    plotter(xx, yy,zz,cpsd,cmsd,c, alpha,beta, case, Z)

def writeOptimalP(case, P, method):
    dirpath = path.join("plots_angle", "optimal_alphabeta.txt")
    #dirpath = path.join("plots_curv", "optimal_alphabeta.txt")
    alpha = P[0]
    beta = P[1]
    #with open(dirpath,"a") as f:
    #    f.write("%s %s alpha=%s beta=%s\n" % (case,method, alpha,beta))


def optimal_angle(i,files,plot):
    files = path.join("param_angle", files)
    
    method = "angle"

    alphabeta_bound = np.loadtxt("alphabeta_bound.txt")


    amin, amax, bmin, bmax = alphabeta_bound[i][0],alphabeta_bound[i][1],alphabeta_bound[i][2],alphabeta_bound[i][3] 

    sd_angle = 19.1

    # Defined SD planes for angle 
    # Tolerance added for adjusting SD 
    # if there are no intersections found
    def apsd(x, y, tol=None):
        if tol is not None:
            return angle0 + sd_angle - tol
        else:
            return angle0 + sd_angle

    def amsd(x, y, tol=None):
        if tol is not None:
            return angle0 - sd_angle + tol 
        else:
            return angle0 - sd_angle  

    def a(x, y, tol=None):
        return angle0

    N = 30
    edge = [amin,amax,bmin,bmax]
    alpha = np.linspace(amin,amax,N)
    beta = np.linspace(bmin,bmax,N) 
    alpha_long = np.linspace(amin,amax,600)
    beta_long = np.linspace(bmin,bmax,600) 
    
    with open(files) as file:
        angle = [[float(digit) for digit in line.split()] for line in file]

    yy, xx = np.meshgrid(beta, alpha)

    points = np.zeros((N,2))
    for i in range(len(xx)):
        points[i] = [alpha[i],beta[i]]

    f = interpolate.interp2d(beta,alpha, angle, kind='cubic')
    angle0 = f(0,0)
    zz = f(beta,alpha)
    
    case = files[-9:-4]
    Z = []
    print "Working on %s" % case
    for plane in [apsd, amsd]: 
        print "Method: %s" % (plane.__name__)
        # Find intersecting points
        # Reduces SD if no points are found
        zeros = intersect(plane,f, alpha_long, beta_long)
        if len(zeros) > 10:
            dx = int(len(zeros) / 10.)
        elif len(zeros) > 0:
            dx = 1
        else:
            empty = True
            tol = 0.1
            maxiter = 50
            it = 0
            print "Found no points..Adjusting SD"
            while empty == True and it < maxiter:
                print "Iterations: %i" % (it + 1)
                zeros = intersect(plane,f, alpha_long, beta_long, tol)
                if len(zeros) > 0:
                    empty = False
                it += 1
                tol += 0.2

            if len(zeros) > 10:
                dx = int(len(zeros) / 10.)
            elif len(zeros) > 0:
                dx = 1

        if len(zeros) > 0:

            points = []
            for p in zeros:
               # if plane.__name__ == "apsd":
               #     if p[1] > 0:
                points.append(p)
                #elif plane.__name__ == "amsd":
                    #if p[1] < 0:
                 #   points.append(p)

            if plane != a:
                P = points[0]
                tol = 0.15
                dist = 10
                for p in points[1:]:
                    dist_tmp = la.norm(np.array(p))
                    if dist_tmp < dist and dist_tmp > tol:
                        dist = dist_tmp
                        P = p
                Z.append(P)
                writeOptimalP(case, P, plane.__name__)

            #plotAlphaBeta(zeros, case, plane.__name__, edge)
            
    
    plus = f(Z[0][1], Z[0][0]) [0]
    minus = f(Z[1][1], Z[1][0]) [0]
    init = angle0[0]
 #   return plus,minus,init

    plotter(xx, yy,zz,apsd,amsd,a, alpha,beta, case, Z)



def plotAlphaBeta(zeros, case, name, edge):
    alpha = np.transpose(zeros)[0]
    beta = np.transpose(zeros)[1]
    
    methods = {"apsd": "Angle + SD", "amsd": "Angle - SD",
            "cpsd": "Curvature + SD", "cmsd": "Curvature - SD"}
    amin = edge[0]
    amax = edge[1]
    bmin = edge[2]
    bmax = edge[3]
    plt.clf()
    #plt.plot(alpha,beta)
    plt.plot(alpha,beta,'ro')
    plt.title("Case %s. Method: %s" % (case, methods[name]))
    plt.xlabel("Alpha")
    plt.ylabel("Beta")
    plt.axis([amin,amax,bmin,bmax])
    plt.xticks(np.arange(amin, amax, 0.05))
    plt.yticks(np.arange(bmin, bmax, 0.05))

    #plt.savefig("./plots_angle/abplot_%s_%s.png" % (case, name), dpi = 250, bbox_inches='tight') 
    #plt.savefig("./plots_curv/abplot_%s_%s.png" % (case, name), dpi = 250, bbox_inches='tight') 
    

def plotter(xx, yy,zz,cpsd,cmsd,c, alpha,beta, case, Z):
    import matplotlib 
    matplotlib.use("TKAGG")
    rc('text', usetex=True)


    label_size = 25

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for t in ax.xaxis.get_major_ticks(): t.label.set_fontsize(20)
    for t in ax.yaxis.get_major_ticks(): t.label.set_fontsize(20)
    for t in ax.zaxis.get_major_ticks(): t.label.set_fontsize(20)
    #rcParams['ztick.labelsize'] = label_size 

    rc('font', family='serif')
    
    sd = 0.045
    #sd = 19.1
    # Plot the surface difference.
    #surf= ax.plot_surface(xx, yy, zz-cpsd(xx,yy),cmap = cm.Spectral, antialiased=False, alpha = 0.6)
    #surf= ax.plot_surface(xx, yy, zz-cmsd(xx,yy),cmap = cm.Spectral, antialiased=False, alpha = 0.6)

    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.set_xlabel(r'$\alpha$',fontsize=25,)
    ax.set_ylabel(r'$\beta$', fontsize=25, )
    #ax.set_zlabel(r'$\theta(\alpha,\beta)$',fontsize=25 )
    ax.set_zlabel(r'$\kappa(\alpha,\beta)$',fontsize=25 )
    ax.xaxis.labelpad=20
    ax.yaxis.labelpad=20
    ax.zaxis.labelpad=20
    ax.tick_params(axis='z', which='major', pad=10)
    ax.tick_params(axis='y', which='major', pad=7)
    ax.tick_params(axis='x', which='major', pad=7)




    # Add a color bar which maps values to colors.
    #fig.colorbar(surf, shrink=0.8, aspect=5)
    p1 = ax.contour(xx, yy, zz-cpsd(alpha,beta) , levels = [0],alpha=0.0)
    p2 = ax.contour(xx, yy, zz-cmsd(alpha,beta) , levels = [0], alpha=0.0)
    p3 = ax.contour(xx, yy, zz-c(alpha,beta) , levels = [0], alpha=0.0)
    try:
        dat1 = p1.allsegs[0][0]
        ax.plot(dat1[:,0], dat1[:,1], sd+0.001 ,color="blue",  linewidth=3, label=r"$\kappa_{Initial} + SD$")
        #ax.plot(dat1[:,0], dat1[:,1], sd+0.001 ,color="blue",  linewidth=3, label=r"$\theta_{Initial} + SD$")
    except:
        pass

    try:
        dat2 = p2.allsegs[0][0]
        #ax.plot(dat2[:,0], dat2[:,1], -sd+0.001, color="red", linewidth=3, label=r"$\theta_{Initial} - SD$")
        ax.plot(dat2[:,0], dat2[:,1], -sd+0.001, color="red", linewidth=3, label=r"$\kappa_{Initial} - SD$")
    except:
        pass
    try:
        dat3 = p3.allsegs[0][0]
        #ax.plot(dat3[:,0], dat3[:,1], 0.001, color="black", linewidth=3, label=r"$\theta_{Initial}$")
        ax.plot(dat3[:,0], dat3[:,1], 0.001, color="black", linewidth=3, label=r"$\kappa_{Initial}$")
    except:
        pass
    #ax.scatter(dat1[:,0], dat1[:,1], sd+0.001, color="k", s=20)
    #ax.scatter(dat2[:,0], dat2[:,1], -sd+0.001, color="k", s=20)
    #ax.scatter(dat3[:,0], dat3[:,1], 0.001, color="k", s=20)
    ax.view_init(azim=100)

    P1 = Z[0]
    P2 = Z[1]
    P3 = [0,0]
    dx = 0.000001
    ax.plot([P1[0],P1[0]+dx], [P1[1],P1[1]+dx], [sd+0.001,sd+0.001+dx], marker="o", linestyle="None", color="black", linewidth=5, label=r"($\alpha, \beta, +SD$)")
    ax.plot([P2[0],P2[0]+dx], [P2[1],P2[1]+dx], [-sd+0.001,-sd+0.001+dx], marker="o", color="black",linestyle="None",  linewidth=5,label=r"($\alpha, \beta, -SD$)")
    ax.plot([P3[0],P3[0]+dx], [P3[1],P3[1]+dx], [0.001,0.001+dx], marker="o", color="red", linestyle="None", linewidth=5, label=r"Origin")
    plt.legend(fontsize=20, ncol=1, fancybox=True,shadow=True, loc="upper right")

    surf= ax.plot_surface(xx, yy, zz-c(xx,yy),cmap = cm.Spectral, alpha = 0.6)
    loc = path.join("plots_angle", "case_%s_angle.png" % (case))
    #loc = path.join("plots_angle", "case_%s_curvature.png" % (case))
    #plt.savefig(loc, bbox_inches='tight')
    plt.show()

def main(plot=True):
    #basedir  = "./param_angle"
    basedir  = "./param_curv"
    files = listdir(basedir)
    files_angle = sorted([f for f in files if f[4] in ["a"] ])
    files_curv = sorted([f for f in files if f[4] in ["c"] ])

    #optimal_angle(0, files_angle[0], plot)
    plus = []
    minus = []
    inits = []
    embed()
    for i in range(len(files_curv)):
        #p,m,init = 
        optimal_curvature(i, files_curv[i], plot)
    """
    for i in range(len(files_angle)):
        p,m,init = optimal_angle(i, files_angle[i], plot)

        plus.append(p)
        minus.append(m)
        inits.append(init)
    cases = sorted(listdir("initialmodels"))
    plus = np.asarray(plus)
    minus = np.asarray(minus)
    inits = np.asarray(inits)
    meanSD(plus,minus,inits, cases)
    """

def meanSD(p, m, init, case):
    N = len(p)

    print "PLUS"
    dS = init - p
    p0 = init
    p1 = p
    mean0 = sum(p0) / len(p0) 
    mean1 = sum(p1) / len(p1) 
    sd0 = np.sqrt(sum([ (xi - mean0)**2 for xi in p0]) / ( len(p0)-1))
    sd1 = np.sqrt(sum([ (xi - mean1)**2 for xi in p1]) / ( len(p1)-1))
    md = sum(dS) / len(dS)
    print "%.2f$\pm$%.2f & %.2f$\pm$%.2f & %.3f\\\\" % ( mean0, sd0, mean1, sd1, md)

    print "MINUS"
    dS = init - m
    p0 = init
    p1 = m
    mean0 = sum(p0) / len(p0) 
    mean1 = sum(p1) / len(p1) 
    sd0 = np.sqrt(sum([ (xi - mean0)**2 for xi in p0]) / ( len(p0)-1))
    sd1 = np.sqrt(sum([ (xi - mean1)**2 for xi in p1]) / ( len(p1)-1))
    md = sum(dS) / len(dS)
    print "%.2f$\pm$%.2f & %.2f$\pm$%.2f & %.3f\\\\" % ( mean0, sd0, mean1, sd1, md)


    for i in range(N):
        ds = init[i] - p[i]

        print "%s & %.3f & %.3f & %.4f\\\\" % (case[i][6:-4], init[i], p[i], ds)

    for i in range(N):
        ds = init[i] - m[i]
        print "%s & %.3f & %.3f & %.4f\\\\" % (case[i][6:-4], init[i], m[i], ds)

main(True)

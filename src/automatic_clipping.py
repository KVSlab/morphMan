from argparse import ArgumentParser
from common import *
import numpy as np
from IPython import embed

def read_command_line():
    """Read arguments from commandline"""
    parser = ArgumentParser()

    parser.add_argument('--d', '--dir_path', type=str, default=".",
                        help="Path to the folder with all the cases")
    parser.add_argument('--m', type=str, default="model_smoothed.vtp", help="Name of the model file")
    parser.add_argument('--c', type=str, default="centerline_complete.vtp", help="Name of the centerline file")
    parser.add_argument('--anu', type=bool, default=False)

    args = parser.parse_args()

    return args.d, args.m, args.c, args.anu


def move_past_sphere(centerline, inlet):
    i = 0 if not inlet else centerline.GetNumberOfPoints() - 1
    j = 0 if inlet else centerline.GetNumberOfPoints() - 1
    r = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(i)
    center = centerline.GetPoint(i)

    MISphere = vtk.vtkSphere()
    MISphere.SetCenter(center)
    MISphere.SetRadius(r*(1./3))
    step = -1 if inlet else 1
    for k in range(i, j, step):
        value = MISphere.FunctionValue(centerline.GetPoint(k))
        value = MISphere.EvaluateFunction(centerline.GetPoint(k))
        print value
        if value >= 0:
            break
    a = centerline.GetPoint(k)
    b = centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(i)

    a = centerline.GetPoint(j)
    b = 0.5*centerline.GetPointData().GetArray(radiusArrayName).GetTuple1(j)
    return a,b


def getBoundingBox(centerline, inlet):
    endPoint = centerline.GetPoint(centerline.GetNumberOfPoints() - 1) if not inlet else centerline.GetPoint(0)

    bottom, bottom_r = move_past_sphere(centerline, inlet)

    centerline = vmtk_centerline_attributes(centerline)#CenterlineAttribiutes(cl)
    centerline = vmtk_centerline_geometry(centerline, smooth=False)#CenterlineAttribiutes(cl)

    E1 = get_array("ParallelTransportNormals", centerline, k=3)
    E1 = E1[E1.shape[0]-1,:]
    T = get_array("FrenetTangent", centerline, k=3)
    T = T[T.shape[0]-1,:]
    E2 = np.zeros(3)

    V = np.eye(3)
    V[:, 0] = T
    V[:, 1] = E1
    V = gram_schmidt(V)

    E1 = V[:,1] * bottom_r * 1.5
    E2 = V[:,2] * bottom_r * 1.6
    T = T * bottom_r * 3 if not inlet else T * bottom_r * 3 * (-1)

    corners = []
    for O in [bottom, bottom + T]:
        for dir1 in [1, -1]:
            for dir2 in [1, -1]:
                corners.append(O + dir1*E1 + dir2*E2)

    #viz(line, [bottom, endPoint] + corners)

    corners = np.array(corners)
    limits = []
    for i in range(3):
        for f in [np.min, np.max]:
            limits.append(f(corners[:,i]))

    return limits


def clipp(dirpath, model_path, centerline_path, anu):
    # Read in initial centerlines 
    centerlines_complete = read_polydata(path.join(dirpath, centerline_path))
    surface = read_polydata(path.join(dirpath, model_path))

    # Create box and clipper objects
    box = vtk.vtkBox()
    clipper = vtk.vtkClipPolyData()
    clipper.SetInputData(surface)
    clipper.SetClipFunction(box)
    
    # Iterate through centerlines 
    ID_list = [0] + range(centerlines_complete.GetNumberOfLines() - anu)
    inlet = True
    for i in ID_list:
        # Find size of clipping box
        # FIXME: Box is placed slightly upstream of outlet
        # FIXME: Nothing at inlet
        limits = getBoundingBox(extract_single_line(centerlines_complete, i), inlet)
        dx = limits[1] - limits[0]
        dy = limits[3] - limits[2]
        dz = limits[5] - limits[4]
        """

        limits[0] += 0.1*dx
        limits[1] -= 0.1*dx
        limits[2] += 0.1*dy
        limits[3] -= 0.1*dy
        limits[4] += 0.1*dz
        limits[5] -= 0.1*dz
        """

        # Set bounds and perform clipping
        box.SetBounds(limits[0], limits[1], limits[2], limits[3], limits[4], limits[5])
        clipper.GenerateClippedOutputOn()
        clipper.InsideOutOn()
        clipper.Update()

        # Map clip to surface model
        filter = vtk.vtkGeometryFilter()
        filter.SetInputData(clipper.GetClippedOutput())
        filter.Update()
        surface.DeepCopy(filter.GetOutput())
        clipper.Update()

        # Rest is outlets
        inlet = False

    write_polydata(surface, "a_test_clipping.vtp")

        #sys.exit(0)

def viz(centerline, points=None):
    """Help method during development to view the results"""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import numpy as np
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    N = centerline.GetNumberOfCells()
    for i in range(N):
        point_ids = vtk.vtkIdList()
        centerline.GetCellPoints(i, point_ids)
        points0 = []
        for k in range(point_ids.GetNumberOfIds()):
            points0.append(centerline.GetPoint(point_ids.GetId(k)))
        arr = np.asarray(points0)
        x = arr[:,0]
        y = arr[:,1]
        z = arr[:,2]
        ax.plot(x, y, z, label=i)
        ax.legend()
        plt.hold("on")
    
    counter = 0
    if points is not None:
        for p in points:
            ax.plot([float(p[0])], [float(p[1])], [float(p[2])], "o", label=N+counter)
            ax.legend()
            plt.hold("on")
            counter += 1
    
    plt.show()


if __name__ == "__main__":
    dir_path, model_path, centerline_path, anu = read_command_line()
    case = "P0134"
    name = "surface"
    dirpath = path.join(dir_path, case, name)
    clipp(dirpath, model_path, centerline_path, anu)

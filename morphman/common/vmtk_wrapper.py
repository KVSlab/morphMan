##   Copyright (c) Aslak W. Bergersen, Henrik A. Kjeldsberg. All rights reserved.
##   See LICENSE file for details.

##      This software is distributed WITHOUT ANY WARRANTY; without even
##      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##      PURPOSE.  See the above copyright notices for more information.

from os import path

from vmtk import vtkvmtk, vmtkscripts

# Global array names
from morphman.common.vtk_wrapper import read_polydata, write_polydata

radiusArrayName = 'MaximumInscribedSphereRadius'
surfaceNormalsArrayName = 'SurfaceNormalArray'
parallelTransportNormalsArrayName = 'ParallelTransportNormals'
groupIDsArrayName = "GroupIds"
abscissasArrayName = 'Abscissas'
branchClippingArrayName = 'BranchClippingArray'


def vmtk_smooth_centerline(centerlines_output, num_iter, smooth_factor):
    centerline_smoothing = vmtkscripts.vmtkCenterlineSmoothing()
    centerline_smoothing.SetInputData(centerlines_output)
    centerline_smoothing.SetNumberOfSmoothingIterations(num_iter)
    centerline_smoothing.SetSmoothingFactor(smooth_factor)
    centerline_smoothing.Update()
    centerlines_output = centerline_smoothing.GetOutput()

    return centerlines_output


def vmtk_compute_centerlines(end_point, inlet, method, outlet, pole_ids, resampling, surface, voronoi):
    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.SeedSelectorName = method
    centerlines.AppendEndPoints = end_point
    centerlines.Resampling = 1
    centerlines.ResamplingStepLength = resampling
    centerlines.SourcePoints = inlet
    centerlines.TargetPoints = outlet
    if voronoi is not None and pole_ids is not None:
        centerlines.VoronoiDiagram = voronoi
        centerlines.PoleIds = pole_ids
    centerlines.Execute()
    centerlines_output = centerlines.Centerlines

    return centerlines, centerlines_output


def vmtk_compute_centerline_sections(surface, centerline):
    """
    Wrapper for vmtk centerline sections.

    Args:
        surface (vtkPolyData): Surface to meassure area.
        centerline (vtkPolyData): centerline to measure along.

    Returns:
        line (vtkPolyData): centerline with the attributes
        centerline_sections_area (vtkPolyData): sections along the centerline
    """
    centerline_sections = vtkvmtk.vtkvmtkPolyDataCenterlineSections()
    centerline_sections.SetInputData(surface)
    centerline_sections.SetCenterlines(centerline)
    centerline_sections.SetCenterlineSectionAreaArrayName('CenterlineSectionArea')
    centerline_sections.SetCenterlineSectionMinSizeArrayName('CenterlineSectionMinSize')
    centerline_sections.SetCenterlineSectionMaxSizeArrayName('CenterlineSectionMaxSize')
    centerline_sections.SetCenterlineSectionShapeArrayName('CenterlineSectionShape')
    centerline_sections.SetCenterlineSectionClosedArrayName('CenterlineSectionClosed')
    centerline_sections.Update()

    centerlines_sections_area = centerline_sections.GetOutput()
    line = centerline_sections.GetCenterlines()

    return line, centerlines_sections_area


def vmtk_compute_geometric_features(line, smooth, outputsmoothed=False, factor=1.0, iterations=100):
    """Wrapper for vmtk centerline geometry.

    Args:
        line (vtkPolyData): Line to compute centerline geometry from.
        smooth (bool): Turn on and off smoothing before computing the geometric features.
        outputsmoothed (bool): Turn on and off the smoothed centerline.
        factor (float): Smoothing factor.
        iterations (int): Number of iterations.

    Returns:
        line (vtkPolyData): Line with geometry.
    """
    geometry = vmtkscripts.vmtkCenterlineGeometry()
    geometry.Centerlines = line

    if smooth:
        geometry.LineSmoothing = 1
        geometry.OutputSmoothedLines = outputsmoothed
        geometry.SmoothingFactor = factor
        geometry.NumberOfSmoothingIterations = iterations
    geometry.FernetTangentArrayName = "FernetTangent"
    geometry.FernetNormalArrayName = "FernetNormal"
    geometry.CurvatureArrayName = "Curvature"
    geometry.TorsionArrayName = "Torsion"
    geometry.Execute()

    return geometry.Centerlines


def vmtk_compute_centerline_attributes(line):
    """ Wrapper for centerline attributes.

    Args:
        line (vtkPolyData): Line to investigate.

    Returns:
        line (vtkPolyData): Line with centerline atributes.
    """
    attributes = vmtkscripts.vmtkCenterlineAttributes()
    attributes.Centerlines = line
    attributes.NormalsArrayName = parallelTransportNormalsArrayName
    attributes.AbscissaArrayName = abscissasArrayName
    attributes.Execute()
    centerlines = attributes.Centerlines

    return centerlines


def vmtk_resample_centerline(line, length):
    """Wrapper for vmtkcenterlineresampling

    Args:
        line (vtkPolyData): line to resample.
        length (float): resampling step.

    Returns:
        line (vtkPolyData): Resampled line.
    """
    resampler = vmtkscripts.vmtkCenterlineResampling()
    resampler.Centerlines = line
    resampler.Length = length
    resampler.Execute()

    line = resampler.Centerlines

    return line


def vmtk_cap_polydata(surface):
    """Wrapper for vmtkCapPolyData

    Args:
        surface (vtkPolyData): Surface to be capped.

    Returns:
        surface (vtkPolyData): Capped surface.
    """
    surface_capper = vtkvmtk.vtkvmtkCapPolyData()
    surface_capper.SetInputData(surface)
    surface_capper.SetDisplacement(0.0)
    surface_capper.SetInPlaneDisplacement(0.0)
    surface_capper.Update()

    return surface_capper.GetOutput()


def vmtk_smooth_surface(surface, method, iterations=800, passband=1.0, relaxation=0.01):
    """Wrapper for a vmtksurfacesmoothing.

    Args:
        surface (vtkPolyData): Input surface to be smoothed.
        method (str): Smoothing method.
        iterations (int): Number of iterations.
        passband (float): The passband for Taubin smoothing.
        relaxation (float): The relaxation for laplace smoothing.

    Returns:
        surface (vtkPolyData): The smoothed surface.

    """
    smoother = vmtkscripts.vmtkSurfaceSmoothing()
    smoother.Surface = surface
    smoother.NumberOfIterations = iterations

    if method == "laplace":
        smoother.RelaxationFactor = relaxation
    elif method == "taubin":
        smoother.PassBand = passband

    smoother.Method = method
    smoother.Execute()
    surface = smoother.Surface

    return surface


def vmtk_compute_voronoi_diagram(surface, filename):
    """
    Creates a surface model's
    coresponding voronoi diagram

    Args:
        surface (vtkPolyData): Surface model
        filename (str): Path where voronoi diagram is stored

    Returns:
        newVoronoi (vtkPolyData): Voronoi diagram
    """
    if path.isfile(filename):
        return read_polydata(filename)

    voronoi = vmtkscripts.vmtkDelaunayVoronoi()
    voronoi.Surface = surface
    voronoi.RemoveSubresolutionTetrahedra = 1
    voronoi.Execute()
    write_polydata(voronoi.VoronoiDiagram, filename)

    new_voronoi = voronoi.VoronoiDiagram

    return new_voronoi


def vmtk_polyball_modeller(complete_voronoi_diagram, poly_ball_size):
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    modeller.SetInputData(complete_voronoi_diagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(poly_ball_size)
    modeller.Update()
    return modeller


def vmtk_surface_connectivity(surface):
    connector = vmtkscripts.vmtkSurfaceConnectivity()
    connector.Surface = surface
    connector.CleanOutput = 1
    connector.Execute()

    return connector


def vmtk_branch_clipper(clipped_centerlines, surface):
    clipper = vmtkscripts.vmtkBranchClipper()
    clipper.Surface = surface
    clipper.Centerlines = clipped_centerlines
    clipper.RadiusArrayName = radiusArrayName
    clipper.GroupIdsArrayName = groupIDsArrayName
    clipper.BlankingArrayName = branchClippingArrayName
    clipper.Execute()

    return clipper


def vmtk_endpoint_extractor(centerlines, clipspheres):
    extractor = vmtkscripts.vmtkEndpointExtractor()
    extractor.Centerlines = centerlines
    extractor.RadiusArrayName = radiusArrayName
    extractor.GroupIdsArrayName = groupIDsArrayName
    extractor.BlankingArrayName = branchClippingArrayName
    extractor.NumberOfEndPointSpheres = clipspheres
    extractor.Execute()

    return extractor


def vmtk_compute_surface_normals(capped_surface):
    surface_normals = vmtkscripts.vmtkSurfaceNormals()
    surface_normals.Surface = capped_surface
    surface_normals.NormalsArrayName = surfaceNormalsArrayName
    surface_normals.Execute()
    capped_surface_with_normals = surface_normals.Surface

    return capped_surface_with_normals


def vmtk_compute_branch_extractor(centerlines_complete):
    brancher = vmtkscripts.vmtkBranchExtractor()
    brancher.Centerlines = centerlines_complete
    brancher.RadiusArrayName = radiusArrayName
    brancher.Execute()
    centerlines_branched = brancher.Centerlines

    return centerlines_branched

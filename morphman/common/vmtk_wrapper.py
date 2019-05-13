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
blankingArrayName = 'Blanking'
branchClippingArrayName = 'BranchClippingArray'


def vmtk_smooth_centerline(centerlines, num_iter, smooth_factor):
    """
    Wrapper for vmtkCenterlineSmoothing. Smooth centerlines with a moving average filter.

    Args:
        centerlines (vtkPolyDat): Centerline to be smoothed.
        num_iter (int): Number of smoothing iterations.
        smooth_factor (float): Smoothing factor

    Returns:
        vtkPolyData: Smoothed version of input centerline
    """
    centerline_smoothing = vmtkscripts.vmtkCenterlineSmoothing()
    centerline_smoothing.SetInputData(centerlines)
    centerline_smoothing.SetNumberOfSmoothingIterations(num_iter)
    centerline_smoothing.SetSmoothingFactor(smooth_factor)
    centerline_smoothing.Update()
    centerlines_smoothed = centerline_smoothing.GetOutput()

    return centerlines_smoothed


def vmtk_compute_centerlines(end_point, inlet, method, outlet, pole_ids, resampling_step, surface, voronoi,
                             flip_normals=False, cap_displacement=None, delaunay_tolerance=None,
                             simplify_voronoi=False):
    """
    Wrapper for vmtkCenterlines.
    compute centerlines from a branching tubular surface. Seed points can be interactively selected on the surface,
    or specified as the barycenters of the open boundaries of the surface.

    Args:
        end_point (int): Toggle append open profile barycenters to centerlines
        surface (vktPolyData): Surface model
        voronoi (vtkPolyData): Voronoi diagram based on previous centerlines (Optional)
        inlet (ndarray): List of source point coordinates
        method (str): Seed point selection method
        outlet (ndarray): List of target point coordinates
        pole_ids (ndarray): Pole ID list of Voronoi diagram (Optional)
        resampling_step (float): Resampling step
        flip_normals (float): Flip normals after outward normal computation
        cap_displacement (float): Displacement of the center points of caps at open profiles along their normals
        delaunay_tolerance (float): Tolerance for evaluating coincident points during Delaunay tessellation
        simplify_voronoi (bool): Toggle simplification of Voronoi diagram

    Returns:

    """
    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.SeedSelectorName = method
    centerlines.AppendEndPoints = end_point
    centerlines.Resampling = 1
    centerlines.ResamplingStepLength = resampling_step
    centerlines.SourcePoints = inlet
    centerlines.TargetPoints = outlet

    if voronoi is not None and pole_ids is not None:
        centerlines.VoronoiDiagram = voronoi
        centerlines.PoleIds = pole_ids
    if flip_normals:
        centerlines.FlipNormals = 1
    if cap_displacement is not None:
        centerlines.CapDisplacement = cap_displacement
    if delaunay_tolerance is not None:
        centerlines.DelaunayTolerance = delaunay_tolerance
    if simplify_voronoi:
        centerlines.SimplifyVoronoi = 1
    centerlines.Execute()
    centerlines_output = centerlines.Centerlines

    return centerlines, centerlines_output


def vmtk_compute_centerline_sections(surface, centerlines):
    """
    Wrapper for vmtk centerline sections.

    Args:
        surface (vtkPolyData): Surface to meassure area.
        centerlines (vtkPolyData): centerline to measure along.

    Returns:
        line (vtkPolyData): centerline with the attributes
        centerline_sections_area (vtkPolyData): sections along the centerline
    """
    centerline_sections = vtkvmtk.vtkvmtkPolyDataCenterlineSections()
    centerline_sections.SetInputData(surface)
    centerline_sections.SetCenterlines(centerlines)
    centerline_sections.SetCenterlineSectionAreaArrayName('CenterlineSectionArea')
    centerline_sections.SetCenterlineSectionMinSizeArrayName('CenterlineSectionMinSize')
    centerline_sections.SetCenterlineSectionMaxSizeArrayName('CenterlineSectionMaxSize')
    centerline_sections.SetCenterlineSectionShapeArrayName('CenterlineSectionShape')
    centerline_sections.SetCenterlineSectionClosedArrayName('CenterlineSectionClosed')
    centerline_sections.Update()

    centerlines_sections_area = centerline_sections.GetOutput()
    line = centerline_sections.GetCenterlines()

    return line, centerlines_sections_area


def vmtk_compute_geometric_features(centerlines, smooth, outputsmoothed=False, factor=1.0, iterations=100):
    """Wrapper for vmtk centerline geometry.

    Args:
        centerlines (vtkPolyData): Line to compute centerline geometry from.
        smooth (bool): Turn on and off smoothing before computing the geometric features.
        outputsmoothed (bool): Turn on and off the smoothed centerline.
        factor (float): Smoothing factor.
        iterations (int): Number of iterations.

    Returns:
        line (vtkPolyData): Line with geometry.
    """
    geometry = vmtkscripts.vmtkCenterlineGeometry()
    geometry.Centerlines = centerlines

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


def vmtk_compute_centerline_attributes(centerlines):
    """ Wrapper for centerline attributes.

    Args:
        centerlines (vtkPolyData): Line to investigate.

    Returns:
        line (vtkPolyData): Line with centerline atributes.
    """
    attributes = vmtkscripts.vmtkCenterlineAttributes()
    attributes.Centerlines = centerlines
    attributes.NormalsArrayName = parallelTransportNormalsArrayName
    attributes.AbscissaArrayName = abscissasArrayName
    attributes.Execute()
    centerlines = attributes.Centerlines

    return centerlines


def vmtk_resample_centerline(centerlines, length):
    """Wrapper for vmtkcenterlineresampling

    Args:
        centerlines (vtkPolyData): line to resample.
        length (float): resampling step.

    Returns:
        line (vtkPolyData): Resampled line.
    """
    resampler = vmtkscripts.vmtkCenterlineResampling()
    resampler.Centerlines = centerlines
    resampler.Length = length
    resampler.Execute()

    resampled_centerline = resampler.Centerlines

    return resampled_centerline


def vmtk_cap_polydata(surface, boundary_ids=None, displacement=0.0, in_plane_displacement=0.0):
    """Wrapper for vmtkCapPolyData.
    Close holes in a surface model.

    Args:
        in_plane_displacement (float): Displacement of boundary baricenters, at section plane relative to the radius
        displacement (float):  Displacement of boundary baricenters along boundary normals relative to the radius.
        boundary_ids (ndarray): Set ids of the boundaries to cap.
        surface (vtkPolyData): Surface to be capped.

    Returns:
        surface (vtkPolyData): Capped surface.
    """
    surface_capper = vtkvmtk.vtkvmtkCapPolyData()
    surface_capper.SetInputData(surface)
    surface_capper.SetDisplacement(displacement)
    surface_capper.SetInPlaneDisplacement(in_plane_displacement)
    if boundary_ids is not None:
        surface_capper.SetBoundaryIds(boundary_ids)
    surface_capper.Update()

    return surface_capper.GetOutput()


def vmtk_smooth_surface(surface, method, iterations=800, passband=1.0, relaxation=0.01, normalize_coordinates=True,
                        smooth_boundary=True):
    """Wrapper for a vmtksurfacesmoothing.

    Args:
        smooth_boundary (bool): Toggle allow change of position of boundary points
        normalize_coordinates (bool): Normalization of coordinates prior to filtering,
            minimize spurious translation effects (Taubin only)
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

    if not normalize_coordinates:
        smoother.NormalizeCoordinates = 0
    if not smooth_boundary:
        smoother.BoundarySmoothing = 0

    smoother.Method = method
    smoother.Execute()
    surface = smoother.Surface

    return surface


def vmtk_compute_voronoi_diagram(surface, filename, simplify_voronoi=False, cap_displacement=None, flip_normals=False,
                                 check_non_manifold=False, delaunay_tolerance=0.001, subresolution_factor=1.0):
    """
    Wrapper for vmtkDelanayVoronoi. Creates a surface model's
    coresponding voronoi diagram.

    Args:
        subresolution_factor (float): Factor for removal of subresolution tetrahedra
        flip_normals (bool): Flip normals after outward normal computation.
        cap_displacement (float): Displacement of the center points of caps at open profiles along their normals
        simplify_voronoi (bool): Use alternative algorith for compute Voronoi diagram, reducing quality, improving speed
        check_non_manifold (bool): Check the surface for non-manifold edges
        delaunay_tolerance (float): Tolerance for evaluating coincident points during Delaunay tessellation
        surface (vtkPolyData): Surface model
        filename (str): Path where voronoi diagram is stored

    Returns:
        new_voronoi (vtkPolyData): Voronoi diagram
    """
    if path.isfile(filename):
        return read_polydata(filename)

    voronoi = vmtkscripts.vmtkDelaunayVoronoi()
    voronoi.Surface = surface
    voronoi.RemoveSubresolutionTetrahedra = 1
    voronoi.DelaunayTolerance = delaunay_tolerance
    voronoi.SubresolutionFactor = subresolution_factor
    if simplify_voronoi:
        voronoi.SimplifyVoronoi = 1
    if cap_displacement is not None:
        voronoi.CapDisplacement = cap_displacement
    if flip_normals:
        voronoi.FlipNormals = 1
    if check_non_manifold:
        voronoi.CheckNonManifold = 1

    voronoi.Execute()
    new_voronoi = voronoi.VoronoiDiagram

    write_polydata(new_voronoi, filename)

    return new_voronoi


def vmtk_polyball_modeller(voronoi_diagram, poly_ball_size):
    """
    Wrapper for vtkvmtkPolyBallModeller.
    Create an image where a polyball or polyball line are evaluated as a function.

    Args:
        voronoi_diagram (vtkPolyData): Input Voronoi diagram representing surface model
        poly_ball_size (list): Resolution of output

    Returns:
        vtkvmtkPolyBallModeller: Image where polyballs have been evaluated over a Voronoi diagram
    """
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    modeller.SetInputData(voronoi_diagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(poly_ball_size)
    modeller.Update()

    return modeller


def vmtk_surface_connectivity(surface, method="largest", clean_output=True, closest_point=None):
    """
    Wrapper for vmtkSurfaceConnectivity. Extract the largest connected region,
    the closest point-connected region or the scalar-connected region from a surface

    Args:
        surface (vtkPolyData): Surface model
        method (str): Connectivity method, either 'largest' or 'closest'
        clean_output (bool): Clean the unused points in the output
        closest_point (ndarray): Coordinates of the closest point

    Returns:
        vmtkSurfaceConnectivity: Filter for extracting largest connected region
    """
    connector = vmtkscripts.vmtkSurfaceConnectivity()
    connector.Surface = surface
    connector.Method = method
    if clean_output:
        connector.CleanOutput = 1
    if closest_point is not None:
        connector.ClosestPoint = closest_point

    connector.Execute()

    return connector


def vmtk_branch_clipper(centerlines, surface, clip_value=0.0, inside_out=False, use_radius_information=True,
                        interactive=False):
    """
    Wrapper for vmtkBranchClipper. Divide a surface in relation to its split and grouped centerlines.

    Args:
        groupIds:
        centerlines (vtkPolyData): Input centerlines
        surface (vtkPolyData): Input surface model
        clip_value (float):
        inside_out (bool): Get the inverse of the branch clipper output.
        use_radius_information (bool): To use MISR info for clipping branches.
        interactive (bool): Use interactive mode, requires user input.

    Returns:
        vmtkBranchClipper: Branch clipper used to divide a surface into regions.
    """
    clipper = vmtkscripts.vmtkBranchClipper()
    clipper.Surface = surface
    clipper.Centerlines = centerlines
    clipper.ClipValue = clip_value
    clipper.RadiusArrayName = radiusArrayName
    clipper.GroupIdsArrayName = groupIDsArrayName
    clipper.BlankingArrayName = blankingArrayName
    if inside_out:
        clipper.InsideOut = 1
    if not use_radius_information:
        clipper.UseRadiusInformation = 0
    if interactive:
        clipper.Interactive = 1

    clipper.Execute()

    return clipper


def vmtk_endpoint_extractor(centerlines, number_of_end_point_spheres, number_of_gap_spheres=1):
    """
    Wrapper for vmtkEndpointExtractor.
    Find the endpoints of a split and grouped centerline

    Args:
        centerlines (vtkPolyData): Input centerlines.
        number_of_end_point_spheres (float): Number of spheres to skip at endpoint
        number_of_gap_spheres (float): Number of spheres to skip per gap.

    Returns:
        vmtkEndpointExtractor: Endpoint extractor based on centerline
    """
    extractor = vmtkscripts.vmtkEndpointExtractor()
    extractor.Centerlines = centerlines
    extractor.RadiusArrayName = radiusArrayName
    extractor.GroupIdsArrayName = groupIDsArrayName
    extractor.BlankingArrayName = branchClippingArrayName
    extractor.NumberOfEndPointSpheres = number_of_end_point_spheres
    extractor.NumberOfGapSpheres = number_of_gap_spheres
    extractor.Execute()

    return extractor


def vmtk_compute_surface_normals(surface, auto_orient_normals=True, orient_normals=True,
                                 compute_cell_normals=False, flip_normals=False):
    """
    Wrapper for vmtkSurfaceNormals.
    Computes the normals of the input surface.

    Args:
        surface (vtkPolyData): Input surface model
        auto_orient_normals (bool): Try to auto orient normals outwards
        orient_normals (bool): Try to orient normals so that neighboring points have similar orientations
        compute_cell_normals (bool): Compute cell normals instead of point normals
        flip_normals (bool): Flip normals after computing them

    Returns:
        vtkPolyData: Surface model with computed normals
    """
    surface_normals = vmtkscripts.vmtkSurfaceNormals()
    surface_normals.Surface = surface
    surface_normals.NormalsArrayName = surfaceNormalsArrayName
    if not auto_orient_normals:
        surface_normals.AutoOrientNormals = 0
    if not orient_normals:
        surface_normals.Consistency = 0
    if compute_cell_normals:
        surface_normals.ComputeCellNormals = 1
    if flip_normals:
        surface_normals.FlipNormals = 1

    surface_normals.Execute()
    surface_with_normals = surface_normals.Surface

    return surface_with_normals


def vmtk_compute_branch_extractor(centerlines):
    """
    Wrapper for vmtkBranchExtractor.
    Split and group centerlines along branches:

    Args:
        centerlines (vtkPolyData): Line to split into branches.

    Returns:
        vtkPolyData: Split centerline.
    """

    brancher = vmtkscripts.vmtkBranchExtractor()
    brancher.Centerlines = centerlines
    brancher.RadiusArrayName = radiusArrayName
    brancher.Execute()
    centerlines_branched = brancher.Centerlines

    return centerlines_branched


def vmtk_surface_curvature(surface, curvature_type="mean", absolute=False,
                           median_filtering=False, curvature_on_boundaries=False,
                           bounded_reciporcal=False, epsilon=1.0, offset=0.0):
    """Wrapper for vmtksurfacecurvature

    Args:
        surface (vtkPolyData): The input surface
        curvature_type (str): The type of surface curvature to compute (mean | gaussian | maximum | minimum)
        absolute (bool): Output the avsolute value of the curvature
        median_filtering (bool): Output curvature after median filtering to suppress numerical noise speckles
        curvature_on_boundaries (bool): Turn on/off curvature on boundaries
        bounded_reciporcal (bool): Output bounded reciprocal of the curvature
        epsilon (float): Bounded reciprocal epsilon at the denominator
        offset (float): Offset curvature by the specified value

    Returns:
        surface (vtkPolydata): Input surface with an point data array with curvature values
    """
    curvature = vmtkscripts.vmtkSurfaceCurvature()
    curvature.Surface = surface
    curvature.CurvatureType = curvature_type
    if absolute:
        curvature.AbsoluteCurvature = 1
    else:
        curvature.AbsoluteCurvature = 0
    if median_filtering:
        curvature.MedianFiltering = 1
    else:
        curvature.MedianFiltering = 0
    if curvature_on_boundaries:
        curvature.CurvatureOnBoundaries = 1
    else:
        curvature.CurvatureOnBoundaries = 0
    if bounded_reciporcal:
        curvature.BoundedReciporcal = 1
    else:
        curvature.BoundedReciporcal = 0
    curvature.Epsilon = epsilon
    curvature.Offset = offset

    curvature.Execute()

    return curvature.Surface

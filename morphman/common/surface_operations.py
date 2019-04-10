def create_new_surface(complete_voronoi_diagram, poly_ball_size=[120, 120, 120]):
    """
    Envelops an input voronoi diagram
    into a new surface model at a
    given resolution determined by
    the poly_ball_size.

    Args:
        complete_voronoi_diagram (vtkPolyData): Voronoi diagram
        poly_ball_size (list): List of dimensional resolution of output model

    Returns:
        envelope (vtkPolyData): Enveloped surface model.
    """
    modeller = vtkvmtk.vtkvmtkPolyBallModeller()
    modeller.SetInputData(complete_voronoi_diagram)
    modeller.SetRadiusArrayName(radiusArrayName)
    modeller.UsePolyBallLineOff()
    modeller.SetSampleDimensions(poly_ball_size)
    modeller.Update()

    # Write the new surface
    marching_cube = vtk.vtkMarchingCubes()
    marching_cube.SetInputData(modeller.GetOutput())
    marching_cube.SetValue(0, 0.0)
    marching_cube.Update()
    envelope = marching_cube.GetOutput()

    return envelope

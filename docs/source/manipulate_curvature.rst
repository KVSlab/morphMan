.. title:: Tutorial: Manipulate curvature

.. _manipulate_curvature:

==============================
Tutorial: Manipulate curvature
==============================

The goal of ``manipulate_curvature.py`` is to increase or decrease the
total curvature/torsion in a vascular segment, see Figure 1 for an example.

.. figure:: Curvature.png
  
  Figure 1: An illustration of the desired output from the method.

In this tutorial, we are using the model with
`ID C0004 <http://ecm2.mathcs.emory.edu/aneuriskdata/download/C0004/C0004_models.tar.gz>`_
from the Aneurisk database. For the commands below we assume that there is a
file `./C0004/surface/model.vtp`, relative to where you execute the command.

In ``manipulate_curvature.py``, there are three options for setting
``region-of-interest``:

 * ``manual``: Manual selection, based on clicking on a surface
 * ``commandline``: Provide the points on the centerline
 * ``first_line``: The section between the inlet and the first bifurcation.

For each point a long the centerline in the region of interest
we need a direction to move the new geometry. To obtain the direction,
we compute the distance between the original and a smoothed
centerline. Using a gaussian smoothing, the new centerline will gradually converge
towards a straight line, depending on the number of iterations and smoothing factor.
We can then easily choose to move each point in the Voronoi diagram correspondingly,
and thus obtain a new surface, like depicted to the left in Figure 2. By simply
inverting the direction of the vector, we can also increase the overall curvature,
see the right most surface in Figure 2.

.. figure:: smoothedsiphon.png

  Figure 2: Sharpened and smoothened version of the siphon.

To reproduce the results shown in Figure 2, run::
        
        python manipulate_curvature.py --ifile C0004/surface/model.vtp --ofile C0004/surface/model.vtp --region-of-interest first_line --smooth_line True

To increase the curvature instead, set ``--smooth_line False``.

For additional information, beyond this tutorial, on the script and
input parameters, please run ``python manipulate_curvature.py -h`` or confer with
the :meth:`curvature_variations`.

.. title:: Miscellaneous

=============
Miscellaneous
=============


.. _landmarking:

Landmarking of the internal carotid artery
==========================================
The internal carotid artery can be classified, or landmarks, into segments.
For an *objective* and *reproduceble* maipulation of the segment
historically referred to as the carotid siphon, we have implemented two
previously published methods for landmarking: Piccinelli et al.
(2011) [1]_, and Bogunović et al. (2014) [2]_ (``--algorithm``).

Although the algorithms are well described in both articles, neither state
how the centerlines are computed, or how the curvature and torsion is derived.
Both algorithms are very sensetive to the input curvature and torision, and
are therefore not directly reproduceble. In Kjeldsberg 2018 [3]_ there is a
a thorough comparison between landmarking algorithms, input parameters,
and centerline smoothing methods which can help you to choose the corret
options for your application.

The script ``automated_landmarking.py`` has three methods for computing
the descrete derivatives of the centerline curve, set with
``--curv-method``.

 1. B-Splines (``spine``)
 2. Discrete derivatives (``disc``)
 3. VMTK (``vmtk``)

To perform landmarking, we will be using the model with `ID C0001 <http://ecm2.mathcs.emory.edu/aneuriskdata/download/C0001/C0001_models.tar.gz>`_
from the Aneurisk database. For the commands below we assume that there is a file `./C0001/surface/model.vtp`
, relative to where you execute the command.
To landmark the surface model, run the following command::

    python automated_landmarking.py --ifile ../test/testdata/C0001/surface/model.vtp --method bogunovic --curv-method spline --curv-method spline --nknots 8

The command will output a file `C0001/surface/landmark_[ALGORITHM]_[CURVMETHOD].particles`
which contains four points defining the interfaces between the segments of the vessel.

.. figure:: landmarking.png

  Figure 1: Landmarked geometry, with interfaces shown as red spheres along the white centerline.


.. _compute_alpha_beta:

Compute alpha and beta
======================
Please see :ref:`manipulate_bend` for a definition of :math:`\alpha` and :math:`\beta`.

Instead of directly setting the extent the model should be moved (``--alpha`` and ``--beta``),
it is more convinient to controll a morphological parameter like maximum curvature, or the
angle in the bend.

The idea behind ``calculate_alpha_beta_values.py`` is to use the centerline as a
proxy for the new geometry, and manipulate only the centerline for a range of ``--alpha`` and
``--beta`` values. The resulting 2D data can be fited to a surface spline, from
which one can easly collect an approporiate value for ``--alpha`` and ``--beta``.

For more information on the input parameters in the script For a more detailed description, please see [3]_.


Common
======
In ``common.py`` we have collected generic functions used by multiple scripts.
Many of the functions wrap existing vtk and vmtk functions in a more pythonic syntax.
Instead of writing 6-7 lines of code to initiate a vtk-object, and set each parameter,
and the input surface, one can call one function with multiple arguments instead,
see for instance :meth:`common.threshold`.

In addition to wrapping vtk and vmtk functionality, there is also new methods for
manipulating centerlines and Voronoi diagrams.

.. [1] Piccinelli, M., Bacigaluppi, S., Boccardi, E., Ene-Iordache, B., Remuzzi, A., Veneziani, A. and Antiga, L., 2011. Geometry of the internal carotid artery and recurrent patterns in location, orientation, and rupture status of lateral aneurysms: an image-based computational study. Neurosurgery, 68(5), pp.1270-1285.
.. [2] Bogunović, H., Pozo, J.M., Cárdenes, R., Villa-Uriol, M.C., Blanc, R., Piotin, M. and Frangi, A.F., 2012. Automated landmarking and geometric characterization of the carotid siphon. Medical image analysis, 16(4), pp.889-903.
.. [3] Kjeldsberg, Henrik Aasen. Investigating the Interaction Between Morphology of the Anterior Bend and Aneurysm Initiation. MS thesis. 2018.

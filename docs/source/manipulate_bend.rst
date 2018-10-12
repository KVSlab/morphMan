.. title:: Tutorial: Manipulate bend

.. _manipulate_bend:

=========================
Tutorial: Manipulate bend
=========================

The goal of ``manipulate_bend.py``, is to alter one specific bend in the
vasculature, like shown in Figure 1.

.. figure:: Bend_moving.png

   Figure 1: Illustration of how we would like to alter a bend.

As shown in Figure 1, we have defined two directions to move a bend in:
:math:`\alpha` and :math:`\beta`. The input parameters ``--alpha`` and
``--beta`` control the distance to move the geometry in each direction,
respectivly. For a detailed description on how we have defined of the
vectors :math:`\alpha` and :math:`\beta`, see [1]_.

In this tutorial, we are using the model with
`ID C0005 <http://ecm2.mathcs.emory.edu/aneuriskdata/download/C0005/C0005_models.tar.gz>`_
from the Aneurisk database. For the commands below we assume that there is a folder
a `./C0005/surface/model.vtp`, relative to where you execute the command.

The region of interest can be set by controlling ``--region-of-interest``
with three different options:

  * ``manual``, provide the start and end of the bend by selecting points interactively
  * ``commandline``, provide two points on the commandline through ``--region-points``
  * ``landmarking``, only valid for the internal carotid artery, execute ``automated_landmarking.py`` prior to running ``manipulate_bend.py``, see :ref:`landmarking`.

Figure 2 depicts an example of modifying the input surface in the :math:`\alpha` ('horizontal') direction only.

.. figure:: bend_alpha_variation.png

  Figure 2: Movement in the vertical direction, determined by :math:`\alpha`.

To recreate the surfaces shown in Figure 2, run the two following commands::

    python manipulate_bend.py --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_vertical_plus.vtp --alpha 0.4  --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 --poly-ball-size 250 250 250

    python manipulate_bend.py --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_vertical_minus.vtp --alpha -0.4  --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 --poly-ball-size 250 250 250

Shown in Figure 3 is the output of changing the surface in the
:math:`\beta` ('vertical') direction only. This can be reproduced by running the two following commands::

    python manipulate_bend.py --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_horizontal_plus.vtp --beta 0.4  --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 --poly-ball-size 250 250 250

    python manipulate_bend.py --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_horizontal_minus.vtp --beta -0.4  --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 --poly-ball-size 250 250 250

.. figure:: bend_beta_variation.png

  Figure 3: Movement in the horizontal direction, determined by :math:`\beta`.

Finally, we can extend the movement to both directions by setting both :math:`\alpha` and :math:`\beta` as command line arguments.
An example, where the bend has been moved in both directions, is illustrated in Figure 4.
In order to reproduce this result, you can run the following commands::

    python manipulate_bend.py --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_plus.vtp --alpha 0.4 --beta 0.4  --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 --poly-ball-size 250 250 250

    python manipulate_bend.py --ifile C0005/surface/model.vtp --ofile C0005/surface/bend_minus.vtp --alpha -0.4 --beta -0.4  --region-of-interest commandline --region-points 49.8 49.7 36.6 53.1 41.8 38.3 --poly-ball-size 250 250 250

.. figure:: bend_alpha_beta_variation.png

  Figure 4: Movement in both directions.

As shown above, the scripts changes the bend as expected, but in general
what would be natural to control is not the distances, but rather a 
morphological parameter like the maximum curvature, or the 'angle' between
the segments. We have therefore added an extra script ``estimate_alpha_and_beta.py``
that can *a priori* find appropriate values for :math:`\alpha` and :math:`\beta` given a
target change in maximum curvature or angle. Please see :ref:`compute_alpha_beta` for more information.

For additional information, beyond this tutorial, on the script and input parameters,
please run ``python manipulate_bend.py -h`` or confer with the :ref:`manipulate_bend`.

.. [1] Kjeldsberg, Henrik Aasen. Investigating the Interaction Between Morphology of the Anterior Bend and Aneurysm Initiation. MS thesis. 2018.

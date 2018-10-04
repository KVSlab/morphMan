.. title:: Tutorial: Manipulate bend

.. _manipulate_bend:

=========================
Tutorial: Manipulate bend
=========================

The goal of ``manipulate_bend.py``, is to alter one specific bend in the vasculature, like shown in Figure 1.

.. figure:: Bend_moving.png

   Figure 1: Illustration of how we would like to alter a bend.

..
    Can be used for a general bend, but if used in ICA...
        Manipulation is initialized by selecting a segment of the vessel, bounded by two clipping points. 
    The two clipping points can be freely chosen along the centerline, but it is highly recommended to landmark the geometry in order to objectively segment the geometry, and use the resulting landmarking points as clipping points.  

Adjusting curvature and angle in the anterior bend utilizes a common script: ``move_siphon.py``. The script performs geometric manipulation of the anterior bend segment, as defined in the landmarking section.
Adjusting the anterior bend relies only on two parameters, the compression/extension factors :math:`\alpha \text{ and } \beta`.
Alteration of the curvature or angle of the anterior bend is performed by specifying these factors in the script  ``automated_geometric_quantities.py`` and ``calculate_alpha_beta_values.py``.
The pipeline for increasing or decreasing either curvature or the bend angle in the anterior bend is described below.   

Alternatively the user may choose any arbitrary values for :math:`\alpha \text{ and } \beta`. 

To perform geometric manipulation of the anterior bend, run the following command::
    
    python move_siphon.py --dir_path [PATH_TO_CASES] --case [CASENAME] --alpha [ALPHA] --beta [BETA]

In general, the compression / extension factors :math:`\alpha \text{ and } \beta` determine the magnitude and direction in which the anterior bend is translated. The manipulation script allows movement in two directions:

* Vertical, determined by :math:`\alpha`
* Horizontal, determined by :math:`\beta`


.. figure:: alpha.png

  Figure 6: Movement in the vertical direction, determined by :math:`\alpha`. FIXME: New model, old model. 


.. figure:: beta.png

  Figure 7: Movement in the horizontal direction, determined by :math:`\beta`. FIXME: New model, old model. 



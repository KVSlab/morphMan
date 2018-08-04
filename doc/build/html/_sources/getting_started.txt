.. title:: Using VMTK Tools

================
Using VMTK Tools
================

Application
===========
This tutorial illustrates the steps for the automated geometric manipulation in patient-spesific morphologies of carotid arteries.  
The framework of the algorithm was originally proposed and implemented in Aslak's and Henrik's master theses (2016, 2018). 

The algorithm relies on the definition of Voronoi Diagram and on its properties, particularly on the fact that given the model surface its Voronoi Diagram can be derived and vice versa.

Initialize cases
================

Prior to manipulation, the models need to be initialized by computing their centerline and Voronoi diagram. 
The goal of the initialization is to achieve three separate centerlines: 

1. Two single centerlines through each daughter vessel
2. A single centerline between the daughter vessels
3. The complete network of centerlines through the model

.. figure:: init_cases.png

  Figure 1: Expected outcome of initialization, using a representative model


To initialize, run the ``initialize_cases.py`` file::

    python initialize_cases.py --dir_path PATH_TO_CASES

Running this commando will cause a render window to pop ut, asking you to specify points on the surface which will act as source points. 
This will trigger the following prompt::

    Please position the mouse and press space to add source points, 'u' to undo

.. figure:: source_pt.png

  Figure 2: Placing source seeds.

To select the source, click on the mesh and a red sphere should appear. 
When the source has been selected, press q to continue. Now you'll be prompted::

   Please position the mouse and press space to add target points, 'u' to undo

.. figure:: target_pt.png

  Figure 3: Placing target seeds.

After selecting target points, press q to continue and the centerlines between the selected points will be computed (may take some time). 
By performing these steps three times, you should have the centerlines specified above.
This concludes the initialization.

Landmarking
===========

Landmarking of the geometry is performed in order to identify different segments of the vessel. 
Identification of specific segments is required in order to perform geometric manipulation. 
The script ``automated_landmarking.py`` includes implementation of two automated landmarking algorithms; the first introduced by Piccinelli et al. (2011), the second by Bogunović et al. (2014). 


Both landmarking algorithms rely on geometric properties of the centerline. 
The script allows the user to select one of four different methods to compute the curvature of the discrete centerline:

1. B-Splines (``spine``)
2. Free-knot regression splines (``freeknot``)
3. Discrete derivatives (``disc``)
4. VMTK (``vmtk``) 

.. todo:: 
    Add references to different splining techniques?

To perform landmarking, run the following command::

    python automated_landmarking.py --dir_path [PATH_TO_CASES] --algorithm [ALGORITHM] --curv_method [METHOD]


This will produce ``landmark_ALGORITHM_CURVMETHOD.particles`` which contains four points defining the interfaces between the segments of the vessel.


.. figure:: landmark.png

  Figure 4: Landmarked geometry, with interfaces shown as red spheres.



Geometric manipulation
======================

The framework presented here allows for geometric manipulation of four independant 
morphological features of the carotid arteries. 

* Bifurcation angle rotation
* Cross-sectional vessel area variation
* Curvature variation in the anterior bend
* Angle variation in the anterior bend

Bifurcation angle rotation
--------------------------
The script ``move_branches.py`` performes an objective rotation of two daughter branches, moving them a given angle along the bifurcation plane. 
In this implementation, positive and negative angles rotate the branches upward and downward, respectively.
The implementation rotates both branches by default, but the user is given the option to only rotate a single branch. 
This is achieved by setting the argument ``leave1`` or ``leave2`` to **True**, to leave the first or second branch, respectively.

To perform rotation of the daughter branches, run the following generalized command::
    
    python move_branches.py --dir_path [PATH_TO_CASES] --case [CASENAME] --angle [ANGLE] --lower True


.. figure:: angle_updown.png

  Figure 5: Rotation of daughter branches, in both directions. 

In addition, the script includes two options for reconstruction of the bifurcation, determined by the parameters ``bif`` and ``lower`` set to **True**.
The ``lower`` parameter creates a more realistic looking bifurcation, while the ``bif`` parameter creates a straight segment between the daughter branches.
A comparison is shown below, where the straight segment is created by running::

    python move_branches.py --dir_path [PATH_TO_CASES] --case [CASENAME] --angle [ANGLE] --bif True

.. figure:: angle_bif.png

  Figure 6: Rotation of daughter branches with different reconstruction of the bifurcation.




Cross-sectional vessel area variation
-------------------------------------
Manipulation of the cross-sectional area of a vessel is performed by running the script ``area_variations.py``. The script performs geometric manipulation along the a vessel, represented by a centerline.
In order to preserve the inlet and the end of the geometry segment, the first and last 10% of the area of interest are left untouched. 

The script ``area_variations.py`` has been extended to include several options related to the vessel area:

* Area variations of the vessel
* Overall area increase / decrease
* Stenosis creation / removal
* Vessel tapering (?) 

Area variations
^^^^^^^^^^^^^^^

Area variation is initialized by a factor :math:`\beta` or a ratio, :math:`R = A_{min} / A_{max}`. 
The script takes ``beta`` and ``ratio`` as command line arguments.
However, the script requires only one of the two parameters to perform area variation. 
Setting :math:`\beta < 0` will cause the ratio :math:`R` to decrease,  whereas :math:`\beta > 0` will cause the ratio to increase. 
The ratio, :math:`R`, behaves as shown in the illustration below. 

To perform area variations of the vessel area, run the following command::
    
    python area_variations.py --dir_path [PATH_TO_CASES] --case [CASENAME] --smooth True --beta 0.5

or::

    python area_variations.py --dir_path [PATH_TO_CASES] --case [CASENAME] --smooth True --ratio 0.5


.. figure:: area_vary.png

  Figure : Area variations throughout the geometry for different ratios. 

Overall area variation
^^^^^^^^^^^^^^^^^^^^^^

The area of interest can also be increase or decreased overall, by using the ``percentage`` argument. 
The ``percentage`` argument determines the percentage to increase or decrease the overall area.

To perform overall increase or decrease of the area of interest, run the following command::
    
    python area_variations.py --dir_path [PATH_TO_CASES] --case [CASENAME] --smooth True --percentage [%]

Below is an illustration of area decrease and increase in a single patient-specific model. 

.. figure:: area_decinc.png

  Figure : Decrease and increase in overall area.

Stenosis creation / removal
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The framework allows for creation or removal of one stenosis located along the area of interest.  
Creation and removal of a stenosis is performed by setting input argument ``stenosis`` to **True**.
Spesifically for stenosis creation, the input arguments ``percentage`` and ``size`` determine the narrowness and length of the stenosis, respectively. 
The ``size`` argument is a factor multiplied with the radius at the selected center point of the stenosis, in order to expand the stenosis exposed area. 
For creation of a stenosis, run the follwing command::
    
    python area_variations.py --dir_path [PATH_TO_CASES] --case [CASENAME] --smooth True --percentage [%] --stenosis True --size [SIZE]

Running this commando will cause a render window to pop ut, asking you to specify points on the surface which will act as the center point of the stenosis. 
This will trigger the following prompt, with a suggested point placement:

.. figure:: single_stenosis.png

  Figure : Placing point where stenosis is centered. 


.. figure:: change_stenosis.png

  Figure : Comparison of new and old model, with and without stenosis.


Similarly, removal of a stenosis is achieved by running the command:: 
    
    python area_variations.py --dir_path [PATH_TO_CASES] --case [CASENAME] --smooth True --stenosis True 

The commando will cause a render window to pop ut, asking you to specify points on the surface which will now act as the boundaries of the stenosis. 
This will trigger the following prompt, with a suggested point placement:

.. figure:: place_stenosis.png

  Figure : Placing points to indentify the boundaries of the stenosis.


.. figure:: fixed_stenosis.png

  Figure : Comparison of new and old model, with and without stenosis.


Curvature and angle variation in the anterior bend
--------------------------------------------------

.. note::
    Manipulation is initialized by selecting a segment of the vessel, bounded by two clipping points. 
    The two clipping points can be freely chosen along the centerline, but it is highly recommended to landmark the geometry in order to objectively segment the geometry, and use the resulting landmarking points as clipping points.  

Adjusting curvature and angle in the anterior bend utilizes a common script: ``move_siphon.py``. The script performs geometric manipulation of the anterior bend segment, as defined in the landmarking section.
Adjusting the anterior bend relies only on two parameters, the compression/extension factors :math:`\alpha \text{ and } \beta`. Selection of these factors in order to increase or decrease the curvature or angle of the anterior bend is performed by the scripts ``automated_geometric_quantities.py`` and ``calculate_alpha_beta_values.py``. The pipeline for increasing or decreasing either peak curvature or the bend angle in the anterior bend is described below.   

Alternatively the user may choose any arbitrary values for :math:`\alpha \text{ and } \beta`. 

To perform geometric manipulation of the anterior bend, run the following command::
    
    python move_siphon.py --dir_path [PATH_TO_CASES] --case [CASENAME] --alpha [ALPHA] --beta [BETA]

In general, the compression / extension factors :math:`\alpha \text{ and } \beta` determine the magnitude and direction in which the anterior bend is translated. The manipulation script allows movement in two directions:

* Vertical, determined by :math:`\alpha`
* Horizontal, determined by :math:`\beta`


.. figure:: alpha.png

  Figure 6: Movement in the vertical direction, determined by :math:`\alpha` 


.. figure:: beta.png

  Figure 7: Movement in the horizontal direction, determined by :math:`\beta` 



Curvature and torsion variation in the vessel
---------------------------------------------

.. figure:: smoothedsiphon.png

  Figure 7: Sharpened and smoothened version of the siphon. 


Selection of compression / extension factors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The compression / extension factors :math:`\alpha \text{ and } \beta` determine the magnitude and direction in which the anterior bend is translated. 
Running the scripts  ``automated_geometric_quantities.py`` and ``calculate_alpha_beta_values.py`` is required if the user is interested in reaching a spesific change in angle or curvature. 





.. title:: Using morphMan

==============
Using morphMan
==============

morphMan (morphological manipulation) is a collection of scripts to objectivly manipulate geometric features
of patient-specific vascular geometries. In this tutorial we examplify the usage
by manipulating `internal carotid arteries <https://en.wikipedia.org/wiki/Internal_carotid_artery/>`_
, but the tool can be applied to any geometry with a tubular shape.

The goal of morphMan is to provide researchers, and other users, with a set of tools to investigate the impact
of altering morphological features in patient-specific geometries, like: area, bifurcation angles, over all
curvature, and a bend. Although these features are possible to manipulate in advaced meshing software,
this is to the best ouf our knowledge, the first objective method for varying a larger set of geometrical
features, while perserving the rest of the geometry.

The geometries used in the tutorials are taken from the Aneurisk repository [1]_,
and are free for anyone to download. You can therefore easly follow the tutorials with the same geometries.

Manipulating surfaces is no trivial task, therefore we derive an alternative representation

All the algorithms are based on lgorithm all utilizes Voronoi Diagrams and on its properties, particularly on the fact that given the model surface its Voronoi Diagram can be derived and vice versa. The algorithms presented in the framework was originally proposed and implemented in Bergersen's and Kjeldsberg's master theses (2016, 2018). TODO: Add link.

.. When the JOSS paper is published, it should be mentioned here
.. When the method paper is published, it should be cited here.

We would like to acknowledge the two open-source projects VTK and VMTK which we base our algorithms on, and to the
authors of Ford et al. 2011 for making the code public, and inspiering us to use manipulate Voronoi diagrams to alter 
surfaces.


Geometric manipulation
======================

The framework presented here allows for geometric manipulation of four independant 
morphological features of the carotid arteries. 

* Bifurcation angle rotation
* Cross-sectional area variation
* Curvature variation in the anterior bend
* Angle variation in the anterior bend

.. [1] AneuriskWeb project website, http://ecm2.mathcs.emory.edu/aneuriskweb. Emory University, Department of Math&CS, 2012.

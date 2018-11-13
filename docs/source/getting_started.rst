.. title:: Using morphMan

.. _getting_started:

==============
Using morphMan
==============

morphMan (morphological manipulation) is a collection of scripts to objectively manipulate morphological features
of patient-specific vascular geometries. In the tutorials, we exemplify the usage
by manipulating `internal carotid arteries <https://en.wikipedia.org/wiki/Internal_carotid_artery>`_
, but morphMan can be applied to any geometry with a tubular shape.

The goal of morphMan is to provide researchers, and other users, with a set of tools to investigate the impact
of altering morphological features in patient-specific geometries. For instance, by combining
morphMan with a computational fluid mechanics solver, it can be used to  *objectively* and
*reproducibly* investigate how cross-sectional area, and overall curvature impacts
the local hemodynamics. Previously, one would have to investigate tens, if not hundreds, of models
to correlate changes in specific morphological parameters to changes in local hemodynamic. Now,
one can keep all other morphological parameters fixed while investigating only the impact of one change.
Thus, morphMan opens a wide range of problems which previously was practically infeasible to study.

Although morphological features are possible to manually manipulate in advanced meshing software,
morphMan is, to the best of our knowledge, the first *objective* and *reproducible* method for
varying morphological features individually, while keeping the remaining geometry unchanged.

The basis for all our algorithms are centerlines and the Voronoi diagram of the surface.
These 'representations' of the surface are easier to manipulate and control than
a surface where all cells are connected. A subset of the algorithms
was presented in Bergersen [1]_ and Kjeldsberg [2]_.

The input surfaces used in the tutorials are taken from the Aneurisk repository [3]_, and are free
for anyone to download. You can therefore easily follow the tutorials with the same geometries as we have used.

We would like to acknowledge the two open-source projects `VTK <https://www.vtk.org>`_
and `VMTK <http://www.vmtk.org>`_, with are the basis of morphMan.


Tutorials
=========
For morphMan to be user friendly, we provide a tutorial for each of the scripts.

* :ref:`manipulate_area` (create/remove a stenosis, in/decrease area variation, in/deflation of a vessel)
* :ref:`manipulate_bend` (Change the curvature or angle of a bend)
* :ref:`manipulate_bifurcation` (Change the angle in a bifurcation)
* :ref:`manipulate_curvature` (Increase or decrease the curvature variation in the a vessel segment)

If you have any questions beyond the tutorials, please do not hesitate to get in touch with us.


New features
============
The existing methods provide many degrees of freedom, however, if you need a specific method
or functionality, please do not hesitate to propose enhancements in the
`issue tracker <https://github.com/KVSlab/morphMan/issues/>`_, or create a pull request with new features.
Our only request is that you follow our
`guidelines for contributing <https://github.com/KVSlab/morphMan/blob/master/CONTRIBUTING.md>`_.

.. [1] Bergersen, Aslak Wigdahl. Investigating the Link Between Patient-specific Morphology and Hemodynamics: Implications for Aneurism Initiation?. MS thesis. 2016.
.. [2] Kjeldsberg, Henrik Aasen. Investigating the Interaction Between Morphology of the Anterior Bend and Aneurysm Initiation. MS thesis. 2018.
.. [3] AneuriskWeb project website, http://ecm2.mathcs.emory.edu/aneuriskweb. Emory University, Department of Math&CS,i 2012.

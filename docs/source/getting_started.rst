.. title:: Using morphMan

==============
Using morphMan
==============

morphMan (morphological manipulation) is a collection of tools to objectivly manipulate morphological features
of patient-specific vascular geometries. In the tutorials we examplify the usage
by manipulating `internal carotid arteries <https://en.wikipedia.org/wiki/Internal_carotid_artery/>`_
, but the tool can be applied to any geometry with a tubular shape.

The goal of morphMan is to provide researchers, and other users, with a set of tools to investigate the impact
of altering morphological features in patient-specific geometries. For instance, morphMan can be used to
investigate how cross-sectional area, bifurcation angles, over all curvature change the local hemodynamics.

Although these features are possible to manually manipulate in advaced meshing software, morphMan is,
to the best ouf our knowledge, the first *objective* and *reproduceble* method for varying a
larger set of morphological features individually, while keeping the rest unchanged.

The basis for all our algorithms are centerlines and the Voronoi diagram of the surface.
These 'representations' of the surface are easier to manipulate and controll than
a surface where all cells are connected. An early version of a subset of the algorithms
were presented in Bergersen [1]_ and Kjeldsberg [2]_.

The geometries used in the tutorials are taken from the Aneurisk repository [3]_, and are free
for anyone to download. You can therefore easly follow the tutorials with the same geometries as we have used.

We would like to acknowledge the two open-source projects `VTK <https://www.vtk.org>`
and `VMTK <http://www.vmtk.org>`, without this tool would not
have been possible. Forthermore, we would also like to acknowledge the authors of Ford et al. 2011 [4]_
for making the code from their publication open-source, and was the starting point for morphMan.


Geometric manipulation
======================

The framework presented enable users to manipulate four geometric features independantly.
Please see the tutorials for additional information:

* :ref:`manipulate_area` (create/remove a stenosis, in/decrease area variation, in/deflation of a vessel)
* :ref:`manipulate_bend` (Change the curvature or angle of a bend)
* :ref:`manipulate_bifurcation` (Change the angle in a bifurcation)
* :ref:`manipulate_curvature` (Increase or decrease the curvature variation in the a vessel segment)


New features
============
These four methods provide many degrees of freedom, however if you need a specific method or functionality, please
do not hesitate to propose enhancements in the `issue tracker <https://github.com/hkjeldsberg/vascularManipulationToolkit/issues/>`_, or create a pull request with new features.

.. [1] Bergersen, Aslak Wigdahl. Investigating the Link Between Patient-specific Morphology and Hemodynamics: Implications for Aneurism Initiation?. MS thesis. 2016.
.. [2] Kjeldsberg, Henrik Aasen. Investigating the Interaction Between Morphology of the Anterior Bend and Aneurysm Initiation. MS thesis. 2018.
.. [3] AneuriskWeb project website, http://ecm2.mathcs.emory.edu/aneuriskweb. Emory University, Department of Math&CS,i 2012.
.. [4] Ford, M.D., Hoi, Y., Piccinelli, M., Antiga, L. and Steinman, D.A., 2009. An objective approach to digital removal of saccular aneurysms: technique and applications. The British Journal of Radiology, 82(special_issue_1), pp.S55-S61.

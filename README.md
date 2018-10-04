morphMan
=====

<p align="center">
    <img src="https://raw.githubusercontent.com/KVSlab/vascularManipulationToolkit/master/docs/source/change_stenosis.png?token=AJKg8GZbIiq2gb9V1elUl0hqm3QrCco4ks5bv20iwA%3D%3D" width="426" height="293" alt="Create a stenois"/>
</p>
<p align="center">
    Manipulation of the Internal Carotid Artery
</p>

Description
-----------

morphMan (morphological manipulation) is a collection of tools to objectivly manipulate morphological features
of patient-specific vascular geometries. In the tutorials we examplify the usage
by manipulating `internal carotid arteries <https://en.wikipedia.org/wiki/Internal_carotid_artery/>`_
, but the tool can be applied to any geometry with a tubular shape.

The goal of morphMan is to provide researchers, and other users, with a set of tools to investigate the impact
of altering morphological features in patient-specific geometries. For instance, morphMan can be used to
investigate how cross-sectional area, bifurcation angles, over all curvature change the local hemodynamics.


Authors
-------

morphMan is developed by

  * Aslak W. Bergersen 
  * Henrik A. Kjeldsberg 

Licence
-------

morphMan is licensed under the GNU GPL, version 3 or (at your option) any
later version.

morphMan is Copyright (2016-2018) by the authors.

Documentation
-------------


If you wish to use morphMan for journal publications, please cite the following master theses: 

* [paper](https://www.duo.uio.no/bitstream/handle/10852/50515/master-bergersen.pdf?sequence=5&isAllowed=y)

* [paper](https://www.duo.uio.no/bitstream/handle/10852/63389/henrikkjeldsberg_master.pdf?sequence=1&isAllowed=y)


moprhMan is specifically designed with complex cylindrical geometries in mind. 


Installation
------------

If you are on Linux or MaxOSX you can install all the packages through anaconda.
First install Anaconda or Miniconda (preferably the python 3.6 version).
Then execute the following command::

       conda create -n your_environment -c vmtk python=3.6 itk vtk vmtk scipy numpy ipython

You can then activate your environment by runing ``source activate your_environment``.
You are now all set and can clone morphMan, and start manipulating your geometries::

    git clone https://github.com/hkjeldsberg/vascularManipulationToolkit


Contact
-------

The latest version of this software can be obtained from

  https://github.com/KVSlab/vascularManipulationToolkit

Please report bugs and other issues through the issue tracker at:
  
  https://github.com/KVSlab/vascularManipulationToolkit/issues

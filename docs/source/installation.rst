.. title:: Installation

============
Installation
============
morphMan is independent scripts which can be run from commandline, and not a python module.
Installation is therefore simply to ensure that the dependencies are installed.
Currently the project is only accessable through `github <https://github.com/hkjeldsberg/vascularManipulationToolkit/>`_.


Compatibility and Dependencies
==============================
The general dependencies of morphMan are 

* VTK > 8.1
* Numpy > 1.13
* SciPy > 1.0.0
* VMTK 1.4

Basic Installation
==================
Please confer with the homepage of each package for installation instructions.
However, if you are on Linux or MaxOSX you can install all the packages through anaconda.
First install Anaconda or Miniconda (preferably the python 3.6 version).
Then execute the following command::

       conda create -n your_environment -c vmtk python=3.6 itk vtk vmtk scipy numpy

You can then activate your environment by runing ``source activate your_environment``.
You are now all set and can clone morphMan, and start manipulating your geometries::

  git clone https://github.com/hkjeldsberg/vascularManipulationToolkit

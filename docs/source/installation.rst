.. title:: Installation

============
Installation
============
morphMan (morphological manipulation) is a collection of tools to objectively manipulate morphological
features of patient-specific vascular geometries. Each tool is a script which can be run from the command line.
Installation is therefore merely to ensure that the dependencies are installed.
Currently, the project is only accessible through `GitHub <https://github.com/hkjeldsberg/vascularManipulationToolkit/>`_.


Compatibility and Dependencies
==============================
The general dependencies of morphMan are 

* VTK > 8.1
* Numpy > 1.13
* SciPy > 1.0.0
* VMTK 1.4

Basic Installation
==================
Please consult with the homepage of each package for installation instructions.
However, if you are on Linux or MaxOSX, you can install all the packages through anaconda.
First, install Anaconda or Miniconda (preferably the python 3.6 version).
Then execute the following command::

  conda create -n your_environment -c vmtk python=3.6 itk vtk vmtk scipy numpy

You can then activate your environment by running ``source activate your_environment``.
You are now all set and can clone morphMan, and start manipulating your geometries::

  git clone https://github.com/KVSlab/morphMan

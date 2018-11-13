.. title:: Installation

============
Installation
============
morphMan (morphological manipulation) is a collection of scripts to objectively manipulate
morphological features of patient-specific vascular geometries. Each script
can be run from the command line, and installation is therefore merely to ensure that all the
dependencies are installed. Currently, the project is only accessible through
`GitHub <https://github.com/KVSlab/morphMan/>`_.


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
However, for simplicty, you can install all the packages (including morphman) through anaconda.
First, install Anaconda or Miniconda (preferably the python 3.6 version).
Then execute the following command::

  conda create -n your_environment -c vmtk python=3.6 itk vtk vmtk scipy numpy morphman

You can then activate your environment by running ``source activate your_environment``.
Now you are all set, and can start using morphMan. morphMan can be accessed by opening a Python console
and typing::

    >>> import morphman

which will give you access to every method within the framework.
Alternatively you can use one of the four main methods of manipulation directly through the terminal, by typing::

    $ morphman-[METHOD]

followed by the command line arguments for the selected method.

Development version
===================
Requirements
~~~~~~~~~~~~
The latest development version of morphMan may be acquired through GitHub.
In order to successfully build and use morphMan, the following additional software needs to be installed on your system:

* Git (>1.6)
* Python (>2.6)

Downloading
~~~~~~~~~~~
The latest development version of morphMan can be found on the official
`morphMan git repository <https://github.com/KVSlab/morphMan>`_ on Github.
Make sure Git (> 1.6) is installed, which is needed to clone the repository.
To clone the morphMan repository, navigate to the directory where you wish
morphMan to be stored, type the following command, and press Enter::

    git clone https://github.com/KVSlab/morphMan

After the source distribution has been downloaded, all the files required will be located
in the newly created ``morphMan`` folder.

Building
~~~~~~~~
In order to build and install morphMan, navigate into the ``morphMan`` folder, where a ``setup.py`` file will be located.
Building and installation of morphMan can be performed by simply running the following command from the terminal window::

    python setup.py install

Usage
~~~~~
You are now ready to use morphMan as a Python module, either in the terminal or in a Python console.
A detailed explanation for usage of morphMan is described in :ref:`getting_started`.

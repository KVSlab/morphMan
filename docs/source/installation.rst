.. title:: Installation

============
Installation
============
morphMan (morphological manipulation) is a collection of scripts to objectively manipulate
morphological features of patient-specific vascular geometries. The project is accessible through
`GitHub <https://github.com/KVSlab/morphMan/>`_ and `Anaconda <https://anaconda.org/morphman/morphman>`_.


Compatibility and Dependencies
==============================
The general dependencies of morphMan are 

* VMTK 1.4.0
* VTK 8.1.0
* Numpy <= 1.13
* SciPy 1.1.0
* Python (2.7 or >=3.5)

Basic Installation
==================
We recommend that you can install morphMan through Anaconda.
First, install Anaconda or Miniconda (preferably the Python 3.6 version).
Then execute the following command in a terminal window::

  $ conda create -n your_environment -c vmtk -c morphman morphman

You can then activate your environment by running ``source activate your_environment``.
Now you are all set, and can start using morphMan. morphMan can be accessed by opening a Python console
and typing::

    >>> import morphman

which will give you access to every method within the framework.
Alternatively you can use one of the six main methods of manipulation directly through the terminal, by typing::

    $ morphman-[METHOD]

followed by the command line arguments for the selected method. A detailed explanation for usage of morphMan is described in :ref:`getting_started`.

.. WARNING:: The VMTK version 1.4, the one currently distributed with Anaconda, has a Python3 bug in `vmtkcenterlines` and `vmtksurfacecurvature`. As a workaround you have to change these files. To find out where it is located please execute::
  
    $ which vmtkcenterlines
    /Users/[Name]/anaconda3/envs/[your_environment]/bin/vmtkcenterlines
    $ python -V
    Python 3.6.2 :: Continuum Analytics, Inc.
  
  Now copy the path up until ``[your_environment]`` and add ``lib/python3.6/site-packages/vmtk/vmtkcenterlines.py``. Please change the path separation symbol to match your operating system and change ``python3.6`` to the python version you are using. If using you are using Miniconda, replace `anaconda3` with `miniconda3`. Using this path you can run the two following lines::

    $ sed -i -e 's/len(self.SourcePoints)\/3/len\(self.SourcePoints\)\/\/3/g' /Users/[Name]/anaconda3/envs/[your_environment]/lib/python3.6/site-packages/vmtk/vmtkcenterlines.py
    $ sed -i -e 's/len(self.TargetPoints)\/3/len\(self.TargetPoints\)\/\/3/g' /Users/[Name]/anaconda3/envs/[your_environment]/lib/python3.6/site-packages/vmtk/vmtkcenterlines.py

  Similarly, for `vmtksurfacecurvature.py`, run the following command::

    $ sed -i -e 's/(len(values) - 1)\/2/\(len\(values\) - 1\)\/\/2/g' /Users/[Name]/anaconda3/envs/[your_environment]/lib/python3.6/site-packages/vmtk/vmtksurfacecurvature.py


Development version
===================

Downloading
~~~~~~~~~~~
The latest development version of morphMan can be found on the official
`morphMan git repository <https://github.com/KVSlab/morphMan>`_ on GitHub.
Make sure Git (>=1.6) is installed, which is needed to clone the repository.
To clone the morphMan repository, navigate to the directory where you wish
morphMan to be stored, type the following command, and press Enter::

    $ git clone https://github.com/KVSlab/morphMan

After the source distribution has been downloaded, all the files required will be located
in the newly created ``morphMan`` folder.

Building
~~~~~~~~
In order to build and install morphMan, navigate into the ``morphMan`` folder, where a ``setup.py``
file will be located. First, make sure that all dependencies are installed. Then, building and installation of morphMan
can be performed by simply running the following command from the terminal window::

    $ python setup.py install


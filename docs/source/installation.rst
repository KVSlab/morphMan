.. title:: Installation
.. highlight:: console

============
Installation
============

morphMan (morphological manipulation) is a collection of scripts to objectively manipulate
morphological features of patient-specific vascular geometries. The project is accessible through
`GitHub <https://github.com/KVSlab/morphMan/>`_ and `Anaconda.org <https://anaconda.org/morphman/morphman>`_.


Compatibility and Dependencies
==============================
The general dependencies of morphMan are 

* VMTK 1.5.0
* VTK 9.1.0
* Numpy 1.23
* SciPy 1.8.1
* Python >=3.8

Basic Installation
==================
We recommend that you can install morphMan through Anaconda.
First, install `Anaconda <https://www.anaconda.com/products/distribution>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_.
Then execute the following command in a terminal window to create a Conda environment with morphMan installed::

  $ conda create -n your_environment -c conda-forge -c morphman morphman

.. note::
    Replace ``your_environment`` with the environment name.

You can then activate your environment by running ``source activate your_environment``.
Now you are all set, and can start using morphMan. morphMan can be accessed by opening a Python console
and typing::

    >>> import morphman

which will give you access to every method within the framework.
Alternatively you can use one of the six main methods of manipulation directly through the terminal, by typing::

    $ morphman-[METHOD]

followed by the command line arguments for the selected method. A detailed explanation for usage of morphMan is described in :ref:`getting_started`.

Development version
===================

Downloading
~~~~~~~~~~~
The latest development version of morphMan can be found on the official
`morphMan git repository <https://github.com/KVSlab/morphMan>`_ on GitHub.
Make sure `Git <https://git-scm.com/>`_ is installed, which is needed to clone the repository.
To clone the morphMan repository, navigate to the directory where you wish
morphMan to be stored, type the following command, and press Enter::

   $ git clone https://github.com/KVSlab/morphMan

If it not already present, this will install Python for you.
After the source distribution has been downloaded, all the files required will be located
in the newly created ``morphMan`` folder.

Building
~~~~~~~~
In order to build and install morphMan, navigate into the ``morphMan`` folder, where a ``setup.py``
file will be located. First, make sure that all dependencies are installed. Then, building and installation of morphMan
can be performed with ``pip`` by running the following command::

    $ python -m pip install .

Alternatively, morphMan can be installed using Python directly (deprecated)::

    $ python setup.py install


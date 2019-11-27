## [morphMan - Morphological Manipulation](https://morphman.readthedocs.io)

[![Documentation Status](https://readthedocs.org/projects/morphman/badge/?version=latest)](https://morphman.readthedocs.io/en/latest/?badge=latest)
[![Build Status](https://travis-ci.org/KVSlab/morphMan.svg?branch=master)](https://travis-ci.org/KVSlab/morphMan)
[![Build status](https://ci.appveyor.com/api/projects/status/2k6q32hqg6g5oopc?svg=true)](https://ci.appveyor.com/project/hkjeldsberg/morphman-s1s38)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01065/status.svg)](https://doi.org/10.21105/joss.01065)

<p align="center">
    <img src="https://raw.githubusercontent.com/KVSlab/morphMan/master/docs/source/make_stenosis.png" width="640 height="280" alt="Create a stenois"/>
</p>
<p align="center">
    Manipulation of the internal carotid artery, by introducing an artificial stenosis.
</p>

Description
-----------
morphMan (morphological manipulation) is a collection of tools for automated and objective 
manipulation and reconstruction of morphological features of patient-specific vascular geometries. 
In the tutorials we examplify the usage
by manipulating [internal carotid arteries](https://en.wikipedia.org/wiki/Internal_carotid_artery)
, but the tool can be applied to any geometry with a tubular shape.

The goal of morphMan is to provide research groups, and other individuals, with a set of tools to investigate the impact
of altering morphological features in patient-specific geometries. For instance, morphMan can be used to
investigate how cross-sectional area, bifurcation angles, and over all curvature change the local hemodynamics.

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
For an introduction to morphMan, and tutorials, please refer to the [documentation](https://morphman.readthedocs.io/en/latest/).

If you wish to use morphMan for journal publications, please cite the journal article from JOSS: 

Kjeldsberg et al., (2019). morphMan: Automated manipulation of vascular geometries. Journal of Open Source Software, 4(35), 1065, https://doi.org/10.21105/joss.01065

Installation
------------

For reference, morphMan requires the following dependencies: VTK > 8.1, Numpy <= 1.13, SciPy > 1.0.0, and VMTK 1.4.
If you are on Windows, macOS or Linux you can install all the general dependencies through anaconda.
First install Anaconda or Miniconda (preferably the Python 3.6 version).
Then execute the following command

        conda create -n your_environment -c vmtk -c morphman morphman

You can then activate your environment by running ``source activate your_environment``.
You are now all set, and can start manipulating your geometries.

Contact
-------
The latest version of this software can be obtained from

  https://github.com/KVSlab/morphMan

Please report bugs and other issues through the issue tracker at:
  
  https://github.com/KVSlab/morphMan/issues

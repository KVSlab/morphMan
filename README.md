## [morphMan - Morphological Manipulation](https://morphman.readthedocs.io)

[![Documentation Status](https://readthedocs.org/projects/morphman/badge/?version=latest)](https://morphman.readthedocs.io/en/latest/?badge=latest)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.01065/status.svg)](https://doi.org/10.21105/joss.01065)
[![codecov](https://codecov.io/gh/KVSlab/morphMan/branch/master/graph/badge.svg)](https://codecov.io/gh/KVSlab/morphMan)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/morphman.svg)](https://anaconda.org/conda-forge/morphman)
[![CI](https://github.com/kvslab/morphman/actions/workflows/lint_and_test.yaml/badge.svg)](https://github.com/kvslab/morphman/actions/workflows/lint_and_test.yaml)


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

morphMan is Copyright (2016-2022) by the authors.

Documentation
-------------
For an introduction to morphMan, and tutorials, please refer to the [documentation](https://morphman.readthedocs.io/en/latest/).

If you use wish to use morphMan for journal publications, please cite the following the following journal article for the methods:

Bergersen, Aslak W., Henrik A. Kjeldsberg, and Kristian Valen‐Sendstad. "A Framework for Automated and Objective Modification of Tubular Structures: Application to the Internal Carotid Artery." International Journal for Numerical Methods in Biomedical Engineering (2020): e3330.

And the journal article from JOSS for the software: 

Kjeldsberg et al., (2019). morphMan: Automated manipulation of vascular geometries. Journal of Open Source Software, 4(35), 1065, https://doi.org/10.21105/joss.01065

Installation
------------

For reference, morphMan requires the following dependencies: VTK 9.1, Numpy 1.23, SciPy 1.8.1, and VMTK 1.5.
If you are on Windows, macOS or Linux you can install all the general dependencies through `conda-forge`.
First install Anaconda or Miniconda, then add `conda-forge` to your channels with:
```
conda config --add channels conda-forge
conda config --set channel_priority strict
```

Once the `conda-forge` channel has been enabled, `morphman` can be installed with `conda`:

```
conda create -n your_environment morphman
```

or with `mamba`:

```
mamba create -n your_environment morphman
```

You can then activate your environment with `morphman` by running 

```
source activate your_environment
```

You are now all set, and can start manipulating your geometries.
To get started, we have created some [tutorials](https://morphman.readthedocs.io/en/latest/getting_started.html) for the different modules within morphMan. 

Contact
-------
The latest version of this software can be obtained from

  https://github.com/KVSlab/morphMan

Please report bugs and other issues through the issue tracker at:
  
  https://github.com/KVSlab/morphMan/issues

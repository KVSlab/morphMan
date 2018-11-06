---
title: 'morphMan: Automated manipulation of vascular geometries'
tags:
  - geometric-algorithms
  - morphology
  - vascular
  - biomedical
authors:
  - name: Henrik A. Kjeldsberg
    affiliation: 1
  - name: Aslak W. Bergersen
    orcid: 0000-0001-5063-3680
    affiliation: 1
  - name: Kristian Valen-Sendstad
    orcid: 0000-0002-2907-0171
    affiliation: 1
affiliations:
 - name: Department of Computational Physiology, Simula Research Laboratory
   index: 1
date: 10 October 2018
bibliography: paper.bib
---

# Summary

Cardiovascular diseases are overwhelming the healthcare systems, and the
costs are anticipated to increase in the years to come [@Murray1997],
not to the mention the personal tragedy for those affected [@gage1996effect].
Systemic risk factors are well known to correlate with cardiovascular diseases in general,
but, for instance, arterial plaques and brain aneurysms are focalized, highlighting
the role of local hemodynamics. Furthermore, blood-flow induced wall shear stress (WSS) is
known to contribute to vessel wall adaption and remodeling [@Malek1999b, @morbiducci2016atherosclerosis],
but is challenging to measure *in-vivo*. On the other hand, medical images are routinely available and have
been extensively used in combination with computational fluid dynamics to
study the initiation, progression, and outcome of vascular pathologies [@taylor2010image].

We know that the morphological features of, for instance, the internal
carotid artery is statistically associated with the presence of aneurysms [@Schimansky01122013, @ingebrigtsen2004bifurcation].
Therefore, understanding how the local hemodynamics change with morphology is of great interest and
is typically investigated with parameterized idealized geometric models [@lauric2018proximal], however at the cost of
oversimplified results. To use realistic geometries we could instead correlate the shape and computed
stresses based on hundreds of patient-specific models, but this is very labor intensive and
error-prone [@valen2018real, @berg2018multiple].

Our goal was to take the best from both approaches and create a tool which could parametrize patient-specific
geometries to mimic the natural variability of morphological features in the population. We here present a framework,
*morphMan*, that allows for *objective*, *reproducible*, and *automatic* virtual manipulation of tubular structures
here exemplified with application to the cerebrovasculature.

![voronoi centerline](./figure1.png)\

**Figure 1:** 
   A visualization of the Voronoi diagram (left) and the centerline (right) of a surface.

In a surface, each cell is connected, and manipulating one will alter the surrounding geometry as well. Instead, we 
have based the algorithms on the centerlines and Voronoi diagram of the surface, see Figure 1. The point in the Voronoi 
diagram are not connected, and therefore easier to manipulate. As a result, only the region of interest is manipulated, and the rest of the geometry is left unchanged. Using the Voronoi diagram to alter the surface was first presented in 
[@Piccinelli2011]; moreover, a subset of the algorithms are presented in [@Bergersen2016] and [@kjeldsberg2018investigating].

In *morphMan v0.1* you can alter cross-sectional area, bifurcation angles, 
overall curvature in a segment, and the shape of specific bends. 
For each
category, there is a wide range of options, thus providing the users with many degrees of
freedom for manipulating the geometries. 
Shown in Figure 2 is an example of rotating the
branches in a bifurcation 'up' and 'down'. 

![manipulation ICA bifurcation](./figure2.png)\

**Figure 2:** 
   Output from *morphMan* for manipulating the bifurcation angles.

Although the framework is applicable by anyone familiar with the Python programming language and basic use of Git/Github, the intended audience of *morphMan* is researchers, particularly within the field of computational biomechanics. Our goal is to provide researchers, and other individuals, with a set of tools to investigate the impact of altering morphological features in patient-specific geometries. For instance, by combining morphMan with a computational fluid mechanics solver, it can be used to objectively and reproducibly investigate how cross-sectional area, and overall curvature impacts the local hemodynamics. Previously, one would have to investigate tens, if not hundreds, of models to correlate changes in specific morphological parameters to changes in local hemodynamic. Now, one can keep all other morphological parameters fixed while investigating only the impact of one change. Thus, morphMan opens a wide range of problems which previously was practically infeasible to study.

To summarize, *morphMan* is a general tool for automated manipulation of tubular objects, easily expandable for specialized manipulations.
In the context of flow in the cardiovascular system,
*morphMan* opens new lines of investigation for unraveling the coupling between
morphology and the computed local hemodynamics.


# Acknowledgements

We acknowledge Alban Souche for testing morphMan, and the two open-source projects [vtk](https://www.vtk.org/), and [vmtk](http://www.vmtk.org).

# References

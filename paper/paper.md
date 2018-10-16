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

Cardiovascular diseases are burdening the healthcare systems and the
costs are expected to rise in the years to come [@Murray1997a].
Systemic risk factors are well known to correlate with cardiovascular diseases in general,
but arterial plaques and brain aneurysms are focalized, highlighting
the role of local hemodynamics. Blood-flow induced wall shear stress (WSS) is
known to contribute to vessel wall adaption and
remodeling [@Malek1999b, @morbiducci2016atherosclerosis], but *in-vivo* measurement of
WSS is challenging. On the other hand, medical images are routinely available and have
been extensively used in combination with computational fluid dynamics to
study the initiation, progression, and outcome of vascular pathologies [@taylor2010image].

We know that the morphological features of, for instance, the internal
carotid artery is statistically associated with the presence of aneurysms [@Schimansky01122013, @ingebrigtsen2004bifurcation].
Therefore, understanding how the local hemodynamics change with morphology is of great interest and
is typically investigated with idealized geometric models [@lauric2018proximal].

The goal of *morphMan* was to develop a framework taking advantage of both
the “patient-specificness” based on medical images and parameterize the natural variability
of morpholigical features to study cerebrovascular flows. We here present a framework that
allows for *objective*, *reproducible*, and *automatic* virtual manipulation of tubular structures,
here exemplified with application to the cerebrovasculature.

The algorithms are based on the centerlines and Voronoi diagram of the surface, see Figure 1. These 'representations'
of the surface are easier to manipulate and control since the cells
are not connected. As a result, the rest of the geometry is left unchanged, and only
the region of interest is manipulated. Using the Voronoi diagram to alter the surface
was first presented in [@Piccinelli2011]; moreover a subset of the algorithms are presented
in [@Bergersen2016] and [@Kjeldsberg2018].

<p align="center">
    <img src="./figure1.png", width="320 height="140" alt="Voronoi diagram and centerline of a model."/>
</p>
<p align="center">
   A visualization of the Voronoi diagram (left) and the centerline (right) of a surface.
</p>

In *morphMan v0.1* you can alter cross-sectional area, bifurcation angles, 
overall curvature in a segment, and the shape of specific bends. For each
category, there is a wide range of options, thus providing the users with many degrees of
freedom for manipulating the geometries. Shown in Figure 2 is an example of rotating
the branches in a bifurcation 'up' and 'down'. *morphMan* is easily expandable for specialized manipulations.

<p align="center">
    <img src="./figure2.png", width="640 height="280" alt="Output of morphMan."\>
</p>
<p align="center">
   Output from *morphMan* for manipulating the bifurcation angles.
</p>

*morphMan* opens new lines of investigation for unraveling the coupling between
morphology and the local hemodynamics, and can guide clinical research and eventually
lead to novel treatment methods.


# Acknowledgements

We acknowledge Alban Souche for testing morphMan, and the two open-source projects [vtk](https://www.vtk.org/), and [vmtk](http://www.vmtk.org).

# References

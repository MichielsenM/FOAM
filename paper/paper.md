---
title: 'Foam: A Python package for Forward Asteroseismic Modelling of gravity modes'
tags:
  - Python
  - astronomy
  - asteroseismology
authors:
  - name: Mathias Michielsen
    orcid: 0000-0001-9097-3655
    # equal-contrib: False
    affiliation: 1 
    # corresponding: true 

affiliations:
 - name: Institute of Astronomy, KU Leuven, Celestijnenlaan 200D, B-3001 Leuven, Belgium
   index: 1
date: 15 September 2023
bibliography: paper.bib
---

# Summary

`Foam` is a python package to perform forward asteroseismic modelling of gravity modes. It automates and streamlines 
a considerable fraction of the modelling process, and can be configured to use various different modelling methodologies.
This includes different ways to match the theoretically predicted oscillations to observations, 
the option to use different sets of observables to compare to their theoretically predicted values, 
the use of different merit functions to determine the goodness of fit, and the option to consider nested subgrids in a statistically meaningfull way.
See @Michielsen2021 and Michielsen et al. (2023) (accepted) for applications of these methodologies to model observed gravity modes.

# Statement of need

Stars spend approximately 90% of their evolution on their so called *main sequence*, during which they fuse hydrogen into helium in their cores. In stars with masses above about 1.2 times the mass of the sun, the stellar core in which these fusion processes take place becomes convective. Macroscopic element transport in and near the
convective cores of these stars has a large influence on their life, since it transports additional hydrogen from outside of the nuclear fusion region into this region. In this way it both prolongs the main-sequence lifetime of stars and enlarges the mass of the helium core at the end of the main sequence, which significantly influences all later stages of their evolution. 
However, these transport processes provide the largest uncertainties in stellar structure and evolution 
models for stars with convective cores, due to our poor understanding of macroscopic element transport and limited number 
of useful observations to test the theories. (See e.g. @Anders2023 for a review on this topic.)

Through asteroseismology, the study of stellar pulsations, we gain the means to unravel
the interior structure of stars [@Aerts2010;@Aerts2021]. Gravity modes in particular have a high sensitivity to 
the properties of the near-core region. We can exploit the probing power of gravity modes, 
observed in e.g. Slowly Pulsating B-type stars [@Waelkens1991], to investigate the physics in the interior of these stars,
particularly the transition region between the convective core and radiative envelope.

`Foam` was developed to streamline the forward modelling process of gravity modes in main-sequence stars with convective cores, 
including multiple options and considerations in the statistical analysis of the grids of stellar equillibrium models.

# Acknowledgements

The research leading to the development of this package has received funding from the Research
Foundation Flanders (FWO) by means of a PhD scholarship to MM under project No. 11F7120N.

# References
---
title: 'Foam: A Python package for forward asteroseismic modelling of gravity modes'
tags:
  - Python
  - astronomy
  - stellar astrophysics  
  - asteroseismology
authors:
  - name: Mathias Michielsen
    orcid: 0000-0001-9097-3655
    affiliation: 1 

affiliations:
 - name: Institute of Astronomy, KU Leuven, Celestijnenlaan 200D, B-3001 Leuven, Belgium
   index: 1
date: 27 September 2023
bibliography: paper.bib
---

# Summary

`Foam` is a python package to perform forward asteroseismic modelling of gravity modes. It automates and streamlines a considerable fraction of the modelling process, comparing grids of stellar models and their oscillation frequencies to sets of frequencies observed in stars. For this purpose, it relies on grids of stellar equilibrium models computed by `MESA` [@Paxton2011;@Paxton2013;@Paxton2015;@Paxton2018;@Paxton2019;@Jermyn2023], and oscillation frequencies computed by `GYRE` [@Townsend2013;@Townsend2018].
`Foam` can be configured to use various different modelling methodologies, including different ways to match the theoretically predicted oscillations to observations, the option to use different sets of observables to compare to their theoretically predicted values, the use of different merit functions to assess the goodness of fit, and the option to consider nested subgrids in a statistically meaningful way. See @Michielsen2021 and @Michielsen2023 for applications of these methodologies to model observed gravity modes. 

# Introduction

Stars spend approximately 90% of their evolution on their so called *main sequence*, during which they fuse hydrogen into helium in their cores. In stars with masses above about 1.2 times the mass of the sun, the stellar core in which these fusion processes take place becomes convective. Macroscopic element transport in and near the convective cores of these stars has a large influence on their life, since it transports additional hydrogen from outside of the nuclear fusion region into this region. In this way it both prolongs the main-sequence lifetime of stars and enlarges the mass of the helium core at the end of the main sequence, which significantly influences all later stages of their evolution. However, these transport processes provide the largest uncertainties in stellar structure and evolution models for stars with convective cores, due to our poor understanding of macroscopic element transport and limited number 
of useful observations to test the theories. [See e.g. @Anders2023 for a review on this topic.]

Through asteroseismology, the study of stellar pulsations, we gain the means to unravel the interior structure of stars [@Aerts2010;@Aerts2021]. Gravity (g-) modes in particular have a high sensitivity to the properties of the near-core region. We can exploit the probing power of g-modes, observed in e.g. Slowly Pulsating B-type stars [@Waelkens1991], to investigate the physics in the interior of these stars, particularly the transition region between the convective core and radiative envelope.

# Statement of need

Some tools have been developed and made publicly available to model and determine stellar parameters of solar-like oscillators, such as `AIMS` [@Rendle2019], `BASTA` [@AguirreBorsenKoch2022], and `pySYD` [@Chontos2022]. However, there are several key differences between the modelling of the pressure (p-) modes observed in solar-like oscillators, and the modelling of the g-modes observed in more massive stars. First and foremost, the well-known asteroseismic scaling relations used for solar-like oscillators cannot be extrapolated to main-sequence stars with a convective core. Secondly, the effect of rotation on p-modes is often included in a perturbative way, whereas the g-mode frequencies are strongly dependent on rotation and require the inclusion of the Coriolis acceleration in a non-perturbative way. Additionally the mass regime of stars with convective cores is subject to strong correlations between several model parameters, which sometimes follow non-linear relationships. To tackle both these non-linear correlations and include uncertainties for the theoretical predictions, the Mahalanobis distance [see @Aerts2018 for its application to asteroseismic modelling] provides a more appropriate merit function than the often used $\chi^2$.

`Foam` was developed to be complimentary to the available modelling tools for solar-like oscillators. It provides a framework for the forward modelling of g-modes in main-sequence stars with convective cores, and tackles the differences in the modelling approach as compared to the case of solar-like oscillators.

# Acknowledgements

The research leading to the development of this package has received funding from the Research
Foundation Flanders (FWO) by means of a PhD scholarship to MM under project No. 11F7120N. MM is grateful to T. Van Reeth for his help concerning the scaling of g-modes with rotation.

# References
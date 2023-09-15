---
layout: homepage
title: Home
permalink: /
---

The main functionality of FOAM is to perform forwards asteroseismic modelling of gravity modes via a modelling pipeline.
An example of a science case for which this pipeline was used is the modelling of the slowly pulsating B-type star KIC7760680 by <a href="https://doi.org/10.1051/0004-6361/202039926" target="_blank"> Michielsen et al. (2021)</a>.
The paper used an early version of the pipeline without optimising the rotation for each model individually, but a lot of its functionality was already used in this paper.
Numerous updates in both performance and functionality have been included since and are employed in the modelling of the B-type binary KIC4930889 by Michielsen et al. (2023) (currently under revision).

There is also some additional functionality outside of the modelling pipeline, mostly in the form of plotting tools, check <a href="https://github.com/MichielsenM/FOAM/blob/master/foam/plot_tools.py" target="_blank"> plot_tools</a> for all available functions. An example to quickly make a Kippenhahn plot or Hertzsprung–Russell diagram from a given <a href="https://docs.mesastar.org/en/latest/index.html" target="_blank"> MESA</a> file in the following way:



<pre><code>import matplotlib.pyplot as plt
from foam import plot_tools

plot_tools.plot_hrd('PATH_TO_MESA_HISTORY')
plot_tools.plot_khd('PATH_TO_MESA_HISTORY')
plt.show()</code></pre>



## Modelling pipeline
The setup is made for the purpose of asteroseismic gravity mode modelling. It's functionality allows for:
- Modelling of the observations with multiple different computed grids.
- Use of different merit functions (reduced chi-squared and Mahalanobis Distance are currently implemented).
- Use of multiple sets of observables. Astroseismic (periods, period spacings, frequencies) but also others (logg, logTeff, logL, ...) can be added in the merit function.
- Use of multiple different methods to construct theoretical pulsation patterns to match the observed pattern.
See [pipeline configuration](./Configuration.md) for a more detailed and complete list of the options and choices.
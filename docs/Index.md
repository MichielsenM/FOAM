---
layout: default
title: Home
permalink: /
---

The main functionality of FoAM is to perform forwards asteroseismic modelling of gravity modes via a modelling pipeline.
There is also some additional functionality outside of the modelling pipeline, mostly in the form of plotting tools, check `foam/plot_tools` for all available functions. An example to quickly make a Kippenhahn plot or Hertzsprung–Russell diagram from a given MESA file in the following way:

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
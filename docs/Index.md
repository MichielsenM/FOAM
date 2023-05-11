---
layout: default
title: Home
permalink: /
---

The main functionality of FoAM is to perform forwards asteroseismic modelling of gravity modes via a [modelling pipeline](#Modelling pipeline).
There is also some [additional functionality](#Functionality outside of the modelling pipeline) outside of the modelling pipeline.


## Modelling pipeline
The setup is made for the purpose of asteroseismic gravity mode modelling. It's functionality allows for:
- Modelling of the observations with multiple different computed grids.
- Use of different merit functions (reduced chi-squared and Mahalanobis Distance are currently implemented).
- Use of multiple sets of observables. Astroseismic (periods, period spacings, frequencies) but also others (logg, logTeff, logL, ...) can be added in the merit function.
- Use of multiple different methods to construct theoretical pulsation patterns to match the observed pattern.
See [pipeline configuration](./Configuration.md) for a more detailed and complete list of the options and choices.

## Functionality outside of the modelling pipeline
### Plotting tools
Some functionality can be used on its own outside of the pipeline, e.g. `foam/plot_tools` could be used to quickly make a Kippenhahn plot or Hertzsprung–Russell diagram from a given MESA file in the following way:

<pre><code>import matplotlib.pyplot as plt
from foam import plot_tools

plot_tools.plot_hrd('PATH_TO_MESA_HISTORY')
plot_tools.plot_khd('PATH_TO_MESA_HISTORY')
plt.show()</code></pre>

### Building grids
The functions in `foam/grid_building_vsc` are an example of how to quickly make a setup for MESA and GYRE grids on the VSC (Vlaams Supercomputer Centrum).

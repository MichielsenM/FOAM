---
layout: default
title: Plotting tools
---

## Plotting tools
There is also some additional plotting functionality outside of the modelling pipeline, you don't need to look at this to do forward asteroseismic modelling, but it can be useful to help analyse MESA output. 
For these plotting tools, check the <a href="https://michielsenm.github.io/FOAM/API/foam/plot_tools.html" target="_blank"> API documentation</a> for all available functions. An example to quickly make a Kippenhahn plot or Hertzsprung–Russell diagram from a given <a href="https://docs.mesastar.org/en/latest/index.html" target="_blank"> MESA</a> file in the following way:


<pre><code>import matplotlib.pyplot as plt
from foam import plot_tools

plot_tools.plot_hrd('PATH_TO_MESA_HISTORY')
plot_tools.plot_khd('PATH_TO_MESA_HISTORY')
plt.show()</code></pre>
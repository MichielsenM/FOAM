# foam

## contents

1. `pipeline`: Pipeline for the asteroseismic modelling of a star, using the functions within the files in this directory.
2. `templates`: A number of template files; GYRE inlists, bash scripts, VSC submit files.
3. `functions_for_gyre.py`: Helpful functions for GYRE input and output.
4. `functions_for_mesa.py`: Helpful functions to process MESA output.
5. `grid_building_slurm.py`: Functions for building MESA and GYRE grids on the SLURM framework.
6. `grid_building_vsc.py`: Functions for building MESA and GYRE grids on the VSC framework.
7. `lambda.csv`: List with lambda (eigenvalue of laplace tidal equations) and nu (spin parameter) for modes up to degree 3. (TAR approximation)
8. `maximum_likelihood_estimator.py`: Functions to perform different kinds of maximum likelihood estimation for the models in a grid, and make correlation plots.
9. `support_functions.py` : Helpful functions for making figures, reading HDF5, processing strings.

### Author
Developed by Mathias Michielsen
```
mathias.michielsen@kuleuven.be
Instituut voor Sterrenkunde
KU Leuven, Belgium
```

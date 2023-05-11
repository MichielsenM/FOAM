# foam

## contents

1. `pipeline`: Pipeline for the asteroseismic modelling of a star, using the functions within the files in this directory.
2. `templates`: A number of template files; GYRE inlists, bash scripts, VSC submit files.
3. `build_optimised_pattern.py`: Selection of the theoretical pulsation patterns that best match the observations, with the possibility to optimise the rotation rate in the process through rescaling
4. `functions_for_gyre.py`: Helpful functions for GYRE input and output.
5. `functions_for_mesa.py`: Helpful functions to process MESA output.
6. `gmode_rotation_scaling`: Calculate g-mode period spacing patterns in the asymptotic regime using the TAR.
7. `grid_building_slurm.py`: Functions for building MESA and GYRE grids on the SLURM framework.
8. `grid_building_vsc.py`: Functions for building MESA and GYRE grids on the VSC framework.
9. `lambda.csv`: List with lambda (eigenvalue of laplace tidal equations) and nu (spin parameter) for modes up to degree 3. (TAR approximation)
10. `maximum_likelihood_estimator.py`: Functions to perform different kinds of maximum likelihood estimation for the models in a grid, and make corner plots.
11. `support_functions.py` : Helpful functions for making figures, reading HDF5, processing strings.

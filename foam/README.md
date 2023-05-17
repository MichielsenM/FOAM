# foam

## contents

1. `pipeline`: Pipeline for the asteroseismic modelling of a star, using the functions within the files in this directory.
2. `additional_constraints`: Extra constraints based on e.g. surface properties of the star, or a binary companion.
3. `build_optimised_pattern.py`: Selection of the theoretical pulsation patterns that best match the observations, with the possibility to optimise the rotation rate in the process through rescaling
4. `functions_for_gyre.py`: Helpful functions for GYRE input and output.
5. `functions_for_mesa.py`: Helpful functions to process MESA output.
6. `gmode_rotation_scaling`: Calculate g-mode period spacing patterns in the asymptotic regime using the TAR.
7. `maximum_likelihood_estimator.py`: Functions to perform different kinds of maximum likelihood estimation for the models in a grid, and make corner plots.
8. `model_grid.py`: Create a summary file to store e.g. all history files of a MESA grid in a nested dictionary.
9. `plot_tools.py`: Plotting functionality
10. `support_functions.py`: Helpful functions for making figures, reading HDF5, processing strings.

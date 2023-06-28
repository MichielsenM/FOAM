---
layout: default
title: Pipeline modules and output
---
# Pipeline modules

The pipeline consists of different modules that are be ran sequentially. Each module produces some intermediate output, enabling modules to be skipped if they have been executed before.
Later steps can thus be repeated with different settings without having to repeat the full pipeline.

By default, most modules will check and skip their execution if the output they would generate is already present. For example output from [grid extraction](#pipe0_extract_grid) is stored in a folder `grid_summary` at the same level as the directories for separate stars, and will only be needed once per (group of) theoretical grid(s). All other output will be stored in the folders for the different modelled stars (see 'Setting up the directory' in [walkthrough](./Walkthrough.md)).
Additionally, when considering nested subgrids, the steps concerning [pattern construction](#pipe1_construct_pattern) need only be performed once for the full non-nested grid without repeating for each nested subgrid.

A lot of the intermediate output is stored in hdf5 files, these can be easily read using the `pandas.read_hdf()` functionality of the [pandas](https://pandas.pydata.org/docs/reference/api/pandas.read_hdf.html) python package.

An example template script for the pipeline is given by `pipeline.py` in `foam/pipeline/`.
The different modules of the pipeline are listed below along with a short explanation of their functionality, and the output they produce. The naming scheme for the output is explained below.

<details>
<summary> <b>Naming scheme clarification</b> (click to expand) </summary> <br>

Each word in the filenames enclosed by {} indicates that it is replaced by a value to indicate certain configuration settings. A short overview is given here:
- `rotation_gyre` the rotational frequency used in the GYRE computations in d^-1
- `kval` meridional degree (k value) of the mode ID (k,m) used in GYRE
- `mval` azimuthal order (m value) of the mode ID (k,m) used in GYRE
- `grid` the name of the specific theoretical grid used.
- `star` the name of the modelled star.
- `observable` indicates if periods or frequencies are used to construct the patterns
- `merit_function` the merit function used to calculate goodness of fit, Mahalanobis Distance is abbreviated to MD, and reduced chi-squared is abbreviated to chi2
- `method` the method used to generate the theoretical frequency patterns to match the observed pattern.
- `suffix_observables` the asteroseismic obsevable used in the merit function. Period, period spacing, and frequency will be abbreviated as P, dP, and f, respectively. Contains '+extra' in case more observables are used in addition to the asteroseismic one. 
- `n_sigma_box` size in standard deviations of the box with models accepted as solutions compatible with the surface properties (logTeff, logL, logg).
</details>

## pipe0_extract_grid
Extract all required information from the MESA profiles and GYRE summary files in the theoretical grids.
Requires the grids to be structured as explained in the [walkthrough](./Walkthrough.md).
This step will be required once for a (group of) theoretical grid(s). Modelling different stars with the same theoretical grid will all use the same summary files that this step creates.

<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

 The `grid_summary/` folder will be created one directory level upwards from the `pipeline.py` script to store
- `surfaceGrid_{grid}.hdf`: the info extracted from the MESA profiles.
- `pulsationGrid_{grid}_rot{rotation_gyre}_k{kval}m{mval}.hdf`: the pulsation information from the GYRE summary files.
</details>

## pipe1_construct_pattern        
Construct the theoretical pulsation patterns, select theoretical pulsation patterns matching the observational pattern, and merge with the models surface properties into one file.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder `extracted_freqs/` to store
- `{observable}_{star}_{grid}_{method}.hdf`: a table with optimised rotation rate (and its error), model parameters, and theoretical frequencies that are matched to the observations.
- `surface+{observable}_{star}_{grid}_{method}.hdf`: the same table combined with the surface properties (logTeff, logL, logg ) of the models.
</details>

After this step, subdirectories will be made in case the pipeline is ran for a nested grid where one or more free parameters are fixed to a certain value.

## pipe2_calculate_likelihood     
Calculate the likelihood of all the theoretical patterns according to the specified merit functions.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder `V_matrix` to store
- `{star}_determinant_conditionNr.tsv`, which holds for each chosen combination of modelling options the condition number of the variance-covariance matrix, and the natural logarithm of the determinant of this matrix (`ln(det(V))`).
- Figures showing the variance-covariance matrices named `{star}_{grid}_{method}_{merit_function}_{suffix_observables}.png`.

Creates folder `meritvalues/` to store
- `{star}_{grid}_{method}_{merit_function}_{suffix_observables}.hdf`: table with the meritvalue assigned by the used merit function, optimised rotation rate, model parameters, and the surface properties (logTeff, logL, logg ...).

</details>

## pipe3_add_constraints
Select all the models that fall inside an n-sigma error box on the provided Teff, logg and logL. If `n_sigma_box` in the configuration is set to None, this step will be skipped. If `constraint_companion` is also provided in the pipeline configuration, those constraints on a binary companion are also taken into account using isochrone-clouds. 
<details>
<summary> <b>Isochrone-clouds</b> (click to expand) </summary> <br>
An isochrone-cloud of a model is made up of all models that have the same metallicity, an age equal within 1 timestep, and whose mass is compatible witin the error margin of the observed mass ratio. (However other parameters can differ between models, e.g. internal mixing processes).
The model of the pulsating star must be compatible with it's observed Teff, logg and logL, while the at least one of the models in its isochrone-cloud must be compatible with the companion's observed Teff, logg and logL.
</details>
<details>

<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder `{n_sigma_box}sigmaBox_meritvalues/` to store
- `{star}_{grid}_{method}_{merit_function}_{suffix_observables}.hdf`: same table as in the [previous step](#pipe2_calculate_likelihood), but only listing the selected models that agree with the n-sigma error box.
(Table with the meritvalue assigned by the used merit function, optimised rotation rate, model parameters, and the surface properties (logTeff, logL, logg ...).)

If constraints of a binary companion are also taken into account, a file `isocloud_grid.h5` will be created holding a nested dictionary.
The nested dictionary will have the grid parameter values as keys, with the order of the nested levels the same as the order of the grid parameters ('free_parameters' followed by 'fixed_parameters', see [pipeline config](./Configuration.md)).
The innermost dictionary will hold certain columns of MESA history files, effectively grouping all the data of a grid that is required to construct isoclouds into one nested dictionary.

</details>

## pipe4_AICc                     
Calculate the [Akaike information criterion (AIC)](https://en.wikipedia.org/wiki/Akaike_information_criterion) corrected for small sample size (AICc).
This AICc is calculated for the best model for each combination of modelling options.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder `{n_sigma_box}sigmaBox_output_tables/` to store
- `{star}_AICc_values_{merit_function}.tsv`: the AICc value of the best model for each chosen combination of modelling options. If the merit function is the Mahalanobis Distance, the condition number of the variance-covariance matrix and the natural logarithm of the determinant of this matrix (`ln(det(V))`) are listed as well.
</details>

## pipe5_best_model_errors        
Calculate the 2 sigma uncertainty region of the maximum likelihood solution using [Bayes' theorem](https://en.wikipedia.org/wiki/Bayes%27_theorem).
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

In folder `{n_sigma_box}sigmaBox_meritvalues/`
- `{star}_{grid}_{method}_{merit_function}_{suffix_observables}_2sigma-error-ellipse.hdf`: same as in [add constraints](#pipe3_add_constraints), but only listing the selected models that fall within the 2 sigma error ellipse according to Bayes' theorem.
</details>

## pipe6_corner_plots             
Make corner plots for all combinations of the different modelling choices.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder `{n_sigma_box}sigmaBox_cornerplots/` to store
- `{star}_{grid}_{method}_{merit_function}_{suffix_observables}.png`: cornerplot with the parameters in the grid and the rotation. The 50% best models are shown, colour-coded according to the log of their merit function value. Models in colour fall within the 2 sigma error ellipse, while those in greyscale fall outside of it. Figures on the diagonal show binned parameter distributions of the models in the error ellipse, and the panel at the top right shows an Hertzsprung-Russell (HR) diagram with 1 and 3 sigma observational error boxes. (The HR diagram is replaced by a Kiel diagram in case the observed logL is not provided but logg is.)
</details>

## pipe7_table_best_models
Write a table with the best model of the grid for each combination of different modelling choices.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

In folder `{n_sigma_box}sigmaBox_output_tables/` 
- `{star}_best-model-table_{merit_function}.txt`: text file containing the best model parameters for each combination of the chosen theoretical grid, seismic observables, and pattern construction methods. These three things are listed first, followed by the grid parameters, the optimal rotation rate of this model, the value of the merit function, and the value of the AICc for this merit function.
</details>

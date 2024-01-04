---
layout: default
title: Pipeline modules and output
---
# Pipeline modules

The pipeline consists of different modules that are ran sequentially. Each module produces some intermediate output, enabling modules to be skipped if they have been executed before.
Later steps can thus be repeated with different settings without having to repeat the full pipeline.

By default, most modules will check and skip their execution if the output they would generate is already present. For example output from [grid extraction](#pipe0_extract_grid) is stored in a folder `grid_summary` at the same level as the directories for separate stars, and will only be needed once per (group of) theoretical grid(s). All other output will be stored in the folders for the different modelled stars (see 'Setting up the directory' in [walkthrough](./Walkthrough.md)).
Additionally, when considering nested subgrids, the steps concerning [pattern construction](#pipe1_construct_pattern) need only be performed once for the full non-nested grid without repeating for each nested subgrid.

A lot of the intermediate output is stored in hdf5 files, these can be easily read using the `pandas.read_hdf()` functionality of the <a href="https://pandas.pydata.org/docs/reference/api/pandas.read_hdf.html" target="_blank"> pandas</a> python package.

An empty template script for the pipeline is given by <a href="https://github.com/MichielsenM/FOAM/tree/master/foam/pipeline/pipeline.py" target="_blank"> foam/pipeline/pipeline.py</a>, and examples with some settings filled in can be found in the
<a href="https://github.com/MichielsenM/FOAM/tree/master/example_setup/KIC7760680/" target="_blank"> example setup</a>.
The different modules of the pipeline are listed below along with a short explanation of their functionality, and the output they produce. The naming scheme for the output is explained below.

<details>
<summary> <b>Naming scheme clarification</b> (click to expand) </summary> <br>

Each word in the filenames enclosed by {} indicates that it is replaced by a value to indicate certain configuration settings. A short overview is given here:
<ul>
<li> <code>rotation_gyre</code> the rotational frequency used in the GYRE computations in d^-1 </li>
<li> <code>kval</code> meridional degree (k value) of the mode ID (k,m) used in GYRE </li>
<li> <code>mval</code> azimuthal order (m value) of the mode ID (k,m) used in GYRE </li>
<li> <code>grid</code> the name of the specific theoretical grid used. </li>
<li> <code>star</code> the name of the modelled star. </li>
<li> <code>observable</code> indicates if periods or frequencies are used to construct the patterns </li>
<li> <code>merit_function</code> the merit function used to calculate goodness of fit, Mahalanobis Distance is abbreviated to MD, and reduced chi-squared is abbreviated to chi2 </li>
<li> <code>method</code> the method used to generate the theoretical frequency patterns to match the observed pattern. </li>
<li> <code>suffix_observables</code> the asteroseismic observable used in the merit function. Period, period spacing, and frequency will be abbreviated as P, dP, and f, respectively. Contains '+extra' in  case more observables are used in addition to the asteroseismic one. </li>
<li> <code>n_sigma_box</code> size in standard deviations of the box with models accepted as solutions compatible with the surface properties (logTeff, logL, logg). </li>
</ul>
</details>

### pipe0_extract_grid
Extract all required information from the MESA profiles and GYRE summary files in the theoretical grids.
Requires the grids to be structured as explained in the [walkthrough](./Walkthrough.md).
This step will be required once for a (group of) theoretical grid(s). Modelling different stars with the same theoretical grid will all use the same summary files that this step creates.

<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

 The <code>grid_summary/</code> folder will be created one directory level upwards from the <code>pipeline.py</code> script to store
 <ul>
<li> <code>surfaceGrid_{grid}.hdf</code>: the info extracted from the MESA profiles. </li>
<li> <code>pulsationGrid_{grid}_rot{rotation_gyre}_k{kval}m{mval}.hdf</code>: the pulsation information from the GYRE summary files. </li>
</ul>
</details>

### pipe1_construct_pattern        
Construct the theoretical pulsation patterns, optimise their rotation rates, select theoretical pulsation patterns matching the observational pattern, and merge with the models surface properties into one file.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder <code>extracted_freqs/</code> to store
<ul>
<li> <code>{observable}_{star}_{grid}_{method}.hdf</code>: a table with optimised rotation rate (and its error), model parameters, and theoretical frequencies that are matched to the observations. </li>
<li> <code>surface+{observable}_{star}_{grid}_{method}.hdf</code>: the same table combined with the surface properties (logTeff, logL, logg ) of the models. </li>
</ul>
</details>

After this step, subdirectories will be made in case the pipeline is ran for a nested grid where one or more free parameters are fixed to a certain value.

### pipe2_calculate_likelihood     
Calculate the likelihood of all the theoretical patterns according to the specified merit functions.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder <code>V_matrix</code> to store
<ul>
<li> <code>{star}_determinant_conditionNr.tsv</code>, which holds for each chosen combination of modelling options the condition number of the variance-covariance matrix, and the natural logarithm of the determinant of this matrix (<code>ln(det(V))</code>). </li>
<li> Figures showing the variance-covariance matrices named <code>{star}_{grid}_{method}_{merit_function}_{suffix_observables}.png</code>. </li>
</ul>
Creates folder <code>meritvalues/</code> to store
<ul>
<li> <code>{star}_{grid}_{method}_{merit_function}_{suffix_observables}.hdf</code>: table with the meritvalue assigned by the used merit function, optimised rotation rate, model parameters, and the surface properties (logTeff, logL, logg ...). </li>
</ul>
</details>

### pipe3_add_constraints
Select all the models that fall inside an n-sigma error box on the provided Teff, logg and logL. If `n_sigma_box` in the configuration is set to None, this step will be skipped. If `constraint_companion` is also provided in the pipeline configuration, those constraints on a binary companion are also taken into account using isochrone-clouds. 
<details>
<summary> <b>Isochrone-clouds</b> (click to expand) </summary> <br>
An isochrone-cloud of a model is made up of all models that have the same metallicity, an age equal within 1 timestep, and whose mass is compatible within the error margin of the observed mass ratio. (However other parameters can differ between models, e.g. internal mixing processes).
The model of the pulsating star must be compatible with it's observed Teff, logg and logL, while at least one of the models in its isochrone-cloud must be compatible with the companion's observed Teff, logg and logL.
</details>
<details>

<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder <code>{n_sigma_box}sigmaBox_meritvalues/</code> to store
<ul>
<li> <code>{star}_{grid}_{method}_{merit_function}_{suffix_observables}.hdf</code>: same table as in the <a href="{{site.baseurl}}/Pipeline#pipe2_calculate_likelihood">previous step</a>, but only listing the selected models that agree with the n-sigma error box.
(Table with the meritvalue assigned by the used merit function, optimised rotation rate, model parameters, and the surface properties (logTeff, logL, logg ...).) </li>
</ul>
If constraints of a binary companion are also taken into account, a file <code>isocloud_grid.h5</code> will be created holding a nested dictionary.
The nested dictionary will have the grid parameter values as keys, with the order of the nested levels the same as the order of the grid parameters (<code>free_parameters</code> followed by <code>fixed_parameters</code>, see <a href="{{site.baseurl}}/Configuration">pipeline config</a>.
The innermost dictionary will hold certain columns of MESA history files, effectively grouping all the data of a grid that is required to construct isoclouds into one nested dictionary.

</details>

### pipe4_AICc
Calculate the <a href="https://en.wikipedia.org/wiki/Akaike_information_criterion" target="_blank"> Akaike information criterion (AIC)</a> corrected for small sample size (AICc).
This AICc is calculated for the best model for each combination of modelling options.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder <code>{n_sigma_box}sigmaBox_output_tables/</code> to store
<ul>
<li> <code>{star}_AICc_values_{merit_function}.tsv</code>: the AICc value of the best model for each chosen combination of modelling options. If the merit function is the Mahalanobis Distance, the condition number of the variance-covariance matrix and the natural logarithm of the determinant of this matrix (<code>ln(det(V))</code>) are listed as well. </li>
</ul>
</details>

### pipe5_best_model_errors        
Calculate the 2 sigma uncertainty region of the maximum likelihood solution using <a href="https://en.wikipedia.org/wiki/Bayes%27_theorem" target="_blank"> Bayes' theorem</a>.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

In folder <code>{n_sigma_box}sigmaBox_meritvalues/</code>
<ul>
<li> <code>{star}_{grid}_{method}_{merit_function}_{suffix_observables}_2sigma-error-ellipse.hdf</code>: same as in <a href="{{site.baseurl}}/Pipeline#pipe3_add_constraints">add constraints</a>, but only listing the selected models that fall within the 2 sigma error ellipse according to Bayes' theorem. </li>
</ul>
</details>

### pipe6_corner_plots             
Make corner plots for all combinations of the different modelling choices.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

Creates folder <code>{n_sigma_box}sigmaBox_cornerplots/</code> to store
<ul>
<li> <code>{star}_{grid}_{method}_{merit_function}_{suffix_observables}.png</code>: cornerplot with the parameters in the grid and the rotation. The 50% best models are shown, colour-coded according to the log of their merit function value. Models in colour fall within the 2 sigma error ellipse, while those in greyscale fall outside of it. Figures on the diagonal show binned parameter distributions of the models in the error ellipse, and the panel at the top right shows an Hertzsprung-Russell (HR) diagram with 1 and 3 sigma observational error boxes. (The HR diagram is replaced by a Kiel diagram in case the observed logL is not provided but logg is.) </li>
</ul>
</details>

### pipe7_table_best_models
Write a table with the best model of the grid for each combination of different modelling choices.
<details>
<summary> <b>Output</b> (click to expand) </summary> <br>

In folder <code>{n_sigma_box}sigmaBox_output_tables/</code>
<ul>
<li> <code>{star}_best-model-table_{merit_function}.txt</code>: text file containing the best model parameters for each combination of the chosen theoretical grid, seismic observables, and pattern construction methods. These three things are listed first, followed by the grid parameters, the optimal rotation rate of this model, the value of the merit function, and the value of the AICc for this merit function. </li>
</ul>
</details>

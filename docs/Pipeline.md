---
layout: default
title: Pipeline modules
---
# Pipeline modules

The pipeline consists of different modules that are be ran sequentially. Each module produces some intermediate output, enabling modules to be skipped if they have been executed before.
Later steps can thus be repeated with different settings without having to repeat the full pipeline.
By default, most modules will check and skip their execution if the output they would generate is already present. For example when considering nested subgrids, the steps concerning [grid extraction](#pipe0_extract_grid) and [pattern construction](#pipe1_construct_pattern) need only be performed once for the full non-nested grid without repeating for each nested subgrid.

A lot of the intermediate output is stored in hdf5 files, these can be easily read using the `pandas.read_hdf()` functionality of the 'pandas' python package.

An example template script for the pipeline is given by `pipeline.py` in `foam/pipeline/`.
The different modules of the pipeline are listed below along with a short explanation of their functionality.

### pipe0_extract_grid
Extract all required information from the MESA profiles and GYRE summary files in the theoretical grids.
Requires the grids to be structured as explained in the [walkthrough](./Walkthrough.md).
It will generate the `grid_summary` folder and store the extracted information there in hdf5 files.

### pipe1_construct_pattern        
Construct the theoretical pulsation patterns, select theoretical pulsation patterns matching the observational pattern, and merge with the models surface properties into one file.

### pipe2_calculate_likelihood     
Calculate the likelihood of all the theoretical patterns according to the specified merit functions.

### pipe3_spectro_constraints      
Select all the models that fall inside an n-sigma spectroscopic error box.

### pipe4_AICc                     
Calculate the [Akaike information criterion (AIC)](https://en.wikipedia.org/wiki/Akaike_information_criterion) corrected for small sample size (AICc).

### pipe5_best_model_errors        
Calculate the 2 sigma uncertainty region of the maximum likelihood solution.

### pipe6_corner_plots             
Make corner plots for all combinations of the different modelling choices.

### pipe7_table_best_models
Write a table with for each different modelling choice the best model of the grid.

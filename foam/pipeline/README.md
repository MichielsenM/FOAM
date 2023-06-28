# Modelling_pipeline
Top level scripts that constitute a modelling pipeline for gravity modes. Run them sequentially as in `pipeline.py`

## contents

1. `pipe0_extract_grid.py`: Extract all required info from MESA profiles and GYRE summary files in a grid.
2. `pipe1_construct_pattern.py`: Construct the theoretical pulsation patterns and merge with surface parameters into one file.
3. `pipe2_calculate_likelihood.py`: Calculate the likelihood of all the theoretical patterns according to the specified merit functions.
4. `pipe3_add_constraints.py`: Select all the models that fall inside an n-sigma error box on surface properties.
5. `pipe4_AICc.py`: Calculate the AICc
6. `pipe5_best_model_errors.py`: Calculate the 2 sigma uncertainty region of the maximum likelihood solution.
7. `pipe6_corner_plots.py`: Make the corner plots of the grid for the different modelling methodologies.
8. `pipe7_table_best_models.py`: Write the best models of the grid as a table

9. `pipeline.py`: The top level script that runs the others sequentially.
10. `pipeline_config.py` : A python class to keep the configuration settings of the modelling pipeline.

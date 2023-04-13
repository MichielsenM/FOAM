# Modelling_pipeline

The setup allows the modelling of a star with the following features:
- Modelling of the observations with multiple computed grids at once.
- With different merit functions (chi-squared and Mahalanobis Distance).
- For multiple sets of observables. Astroseismic (periods, period spacings, frequencies) but also others (logg, logTeff, logL) can go in the merit function.
- For multiple methods to construct the theoretical pulsation patterns.

## How to use this pipeline
Make a folder where you want to have all the output files, and copy `pipeline.py` to it. This script will run all the others, but can also be modified to rerun later parts when the earlier ones have been run at least once.
Copy `example_config.py` to the folder from where you run the `pipeline.py` script and rename to 'config.py'. Adjust this file to configure the different aspects of the pipeline.
Add a file with the observations in a .tsv format. Its columns are the observables and their errors (suffix '\_err' for the error corresponding to a column). Frequencies, periods, and their errors are required for the full functionality of the pipeline. Note that effective temperature should be added as Teff, but will be processed as logTeff (since L and g are normally used in log as well).
Frequencies are preferred to be listed decreasing in value, so increasing in period. The scripts can handle the other way around if the theoretical grid is also extracted in that way. (If the frequencies are ordered decreasing in value, the highest frequency parts in an interrupted pattern will be matched first, which is preferred since those have a lower relative uncertainty). Missing frequencies in the observed pattern can be indicated by 'f_missing' as index and 0 values for the frequency, period, and their errors.

## contents

1. `pipe0_extract_puls_and_spectro.py`: Extract all spectroscopic info from MESA profiles and all frequencies from GYRE models in a grid.
2. `pipe1_construct_pattern.py`: Construct the theoretical pulsation patterns and merge with spectroscopic info into one file.
3. `pipe2_calculate_likelihood.py`: Calculate the likelihood of all the theoretical patterns according to the specified merit functions.
4. `pipe3_spectroConstraints.py`: Select all the models that fall inside an n-sigma spectroscopic error box.
5. `pipe4_AICc.py`: Calculate the AICc
6. `pipe5_bestModel_errors.py`: Calculate the 2 sigma uncertainty region of the maximum likelihood solution.
7. `pipe6_correlationPlots.py`: Make the correlation plots of the grid for the different modelling methodologies.
8. `pipe7_table_bestModels.py`: Write the best models of the grid as a table

9. `pipeline.py`: The top level script that runs the others sequentially.
10 `pipelineConfig.py` : A python class to keep the configuration settings of the modelling pipeline.

### Author
Developed by Mathias Michielsen
```
mathias.michielsen@kuleuven.be
Instituut voor Sterrenkunde
KU Leuven, Belgium
```

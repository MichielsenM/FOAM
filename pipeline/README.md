# Modelling_pipeline

The setup allows the modelling of a star with the following features:
- Modelling of the observations with multiple computed grids at once.
- With different merit functions (chi-squared and Mahalanobis Distance).
- For multiple sets of observables. Astroseismic (periods, period spacings, those 2 combined, frequencies) but also others (logg, logTeff, logL) can go in the merit function.
- For multiple methods to construct the theoretical pulsation patterns.

Make a folder where you want to have all the output files, and copy `pipeline.py` to it. This script will run all the others, but can also be modified to rerun later parts when the earlier ones have been run at least once.
Copy `example_config.py` to the folder from where you run the `pipeline.py` script and rename to 'config.py'. Adjust this file to configure the different aspects of the pipeline.
Add a file with the observations in a .tsv format. Its columns are the observables and their errors (suffix '\_err' for the error corresponding to a column). Note that effective temperature should be added as Teff, but will be processed as logTeff (since L and g are normally used in log as well).
Frequencies should be increasing in value (so decreasing in period). Missing frequencies in the pattern can be indicated by 'f_missing' as index and 0 values for the frequency, period, and their errors.

## contents

1. `0_extract_puls&spectro.py`: Extract all spectroscopic info from MESA profiles and all frequencies from GYRE models in a grid.
2. `1_constuct_pattern.py`: Construct the theoretical pulsation patterns and merge with spectroscopic info into one file.
3. `2_calculate_likelihood.py`: Calculate the likelihood of all the theoretial patterns according to the specified merit functions.
4. `3_spectroClip_AICc.py`: Calculate the AICc, ignoring all the models that fall outside an n-sigma spectroscopic error box.
5. `4_bestModel_errors.py`: Calculate the 2 sigma uncertainty region of the maximum likelihood solution.
6. `5_correlationPlots.py`: Make the correlation plots of the grid for the different modelling methodologies.

7. `example_config.py`: Example of a configuration file for the pipeline
8. `pipeline.py`: The top level script that runs the others sequentially.
9. `table_bestModels.py`: Write the best models of the grid as a LaTeX table.

### Author
Developed by Mathias Michielsen
```
mathias.michielsen@kuleuven.be
Instituut voor Sterrenkunde
KU Leuven, Belgium
```

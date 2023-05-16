---
layout: default
title: Configuration
---
# Pipeline configuration

The modelling pipeline can be configured for various options and modelling choices.
These are set via keyword arguments, supplied upon initialisation of the instance of the PipelineConfig() class.
In the example pipeline script `foam/pipeline/pipeline.py`, this instance is initialised in the beginning of the script.
The following lines are copied from the pipeline, but some of the keyword arguments have already been filled in:
<pre><code>from foam.pipeline import pipeline_config
pipeline_config.config = pipeline_config.PipelineConfig(star='KIC000',
                                                        observations='data_KIC000.tsv',
                                                        rotation_gyre = 0.5)
</code></pre>

All keyword arguments are listed below, grouped in categories.

### Observational data
- star
>   type: string <br>
    default: None <br>
    Name of the star, used for generating filenames

- observations
>   type: string <br>
    default: None <br>
    Full path to the file with the observational data

- periods_or_frequencies_observed
>   type: string <br>
    default: 'period' <br>
    options: 'period' or 'frequency' <br>
    Use the observed periods or frequencies, be consistent in observable_list later.

- highest_amplitude_pulsation
>   type: dictionary of list <br>
    default: None <br>
    Only required when 'highest_amplitude' method is included in 'pattern_methods' <br>
    Pulsation with the highest amplitude to build pattern from when using 'highest_amplitude' method. <br>
    List with highest amplitudes per part of the split pattern,
    ordered the same as the file with the observations. (List of lenght 1 in case of continuous pattern.)

### Simulated theoretical model grid
- grid_parent_directory
>   type: string <br>
    default: None <br>
    Parent directory of the computed grids

- grids
>   type: list of strings <br>
    default: None <br>
    Names of the subdirectories of the grids with different physics

- subgrid
>   type: string <br>
    default: '*' <br>
    String to only select models that have this string in their name. Mainly for testing purposes.
    This is used to create the grid summary files for a small subgrid of the original grid,
    e.g. for making a smaller grid with less free parameters to run some quick tests.

### GYRE and mode information
- kval
>   type: int <br>
    default: 0 <br>
    meridional degree (k value) of the mode ID (k,m) of the g-mode pattern

- mval
>   type: int <br>
    default: 1 <br>
    azimuthal order (m value) of the mode ID (k,m) of the g-mode pattern

- rotation_gyre
>   type: float <br>
    default: 0.0 <br>
    Rotation rate in d^-1 (CYC_PER_DAY of GYRE) that was used in GYRE to calculate the pulsation frequencies

- gyre_dir
>   type: string <br>
    default: os environment variable $GYRE_DIR <br>
    Path to the GYRE directory.

### Modelling methodology
- pattern_methods
>   type: list of strings <br>
    options: 'highest_amplitude', 'highest_frequency', 'chisq_longest_sequence' <br>
    default: ['highest_amplitude', 'highest_frequency', 'chisq_longest_sequence'] <br>
    List of methods to construct the theoretical frequency pattern (repeats modelling for each method). <br>
    >> 'highest_amplitude' builds the pattern starting from the mode with the highest observed amplitude.
    >
    >> 'highest_frequency' will begin matching mode periods starting from the theoretical
    period that is closest to the highest-frequency mode detected in the observed pattern.
    >
    >>  'chisq_longest_sequence' will match each observed mode period to its best matching
        theoretical counterpart, and adopt the longest sequence of consecutive modes retrieved in this way.
        In the case of multiple mode series with the same length, a final pattern
        selection is made based on the best match between theory and observations.

- merit_functions
>   type: list of strings <br>
    options: 'chi2', 'mahalanobis' <br>
    default: ['chi2', 'mahalanobis'] <br>
    List of merit functions (repeats modelling for each merit function).
    >> 'chi2' uses a reduced chi-squared function
    >
    >> 'mahalanobis' uses Mahalanobis distances, see [Aerts et al. 2018](https://ui.adsabs.harvard.edu/abs/2018ApJS..237...15A/abstract)
    for its application in the context of asteroseismology.

- observable_list
>   type: list of list of strings <br>
    default: [['period'], ['period_spacing']] <br>
    Lists of observables to fit (repeats modelling for each list) <br>
    e.g. ['period'] can be expanded to ['period', 'logg'] to include more observables in the likelihood estimation
    Needs to include 'period' or 'period_spacing' if 'period' is used in keyword
    'periods_or_frequencies_observed', and needs to include 'frequency' if 'frequency' is used in that keyword.

- observable_aic
>   type: list of strings <br>
    default: ['P', 'dP'] <br>
    calculate AICc for these observables (abbreviated names) #TODO

- n_sigma_spectrobox
>   type: int <br>
    default: 3 <br>
    Ignore models outside of the n-sigma spectroscopic error box, set to None to include all models.

- free_parameters
>   type: list <br>
    default: ['Z', 'M', 'logD', 'aov', 'fov', 'Xc'] <br>
    List of varied parameters in the grid, as written in grid filenames, (see [Walkthrough](./Walkthrough))
    that remain free parameters in the modelling.

- fixed_parameters
>   type: dict <br>
    default: None <br>
    Dictionary with varied parameters in the grid (dict keys) that are fixed to a certain value (dict value)
    to consider nested grids with fewer free parameters. Defaults to None to not fix any parameters.
    (E.g. {'aov': 0.0, 'fov': 0.0} to fix both these mixing parameters to zero.)

- N_periods
>   type: int <br>
    default: None <br>
    Number of periods in the observed pattern

- N_pattern_parts
>   type: int <br>
    default: 1 <br>
    In how many parts the observed pattern is split. Defaults to 1 assuming an uninterrupted pattern.

### Modelling binaries and enforcing constraints of the companion star
- spectro_companion
>   type: dictionary <br>
    default: None <br>
    Defaults to None instead of a dictionary to not include any constraints from binarity.
    Dictionary with the following keys holding all binary Information (q = Mass_secondary/Mass_primary) e.g.  <br> {'q': \<float\>, 'q_err': \<float\>,'Teff': \<float\>, 'Teff_err': \<float\>,
          'logg': \<float\>, 'logg_err': \<float\>, 'logL': \<float\>, 'logL_err':\<float\>, 'primary_pulsates':\<boolean\>}  <br>
    Set one of the observables Teff, logg, logL to None to not use that observable.


- isocloud_grid_directory
>   type: string <br>
    default: None <br>
    The path to the isocloud grid directory.


### Other settings
- debugging
>   type: boolean <br>
    default: False <br>
    Set to True to set logger level to debug    

- nr_cpu
>   type:int <br>
    default: None <br>
    Number of worker processes to use in multiprocessing. The default 'None' will cause the pools to use the number returned by os.cpu_count().

- conerplot_axis_labels
>   type: dictionary, keys and values are strings <br>
    default: ` {'rot': r'$\Omega_{\mathrm{rot}}$ [d$^{-1}$]' ,'M': r'M$_{\rm ini}$', 'Z': r'Z$_{\rm ini}$', 'logD':r'log(D$_{\rm env}$)', 'aov':r'$\alpha_{\rm CBM}$','fov':r'f$_{\rm CBM}$','Xc':r'$\rm X_c$'} ` <br>
    keys are the grid parameters, values are how they should be put in the labels of the cornerplots' axis

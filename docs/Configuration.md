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
    Name of the star, used for generating filenames. Do not put underscores (`_`) in the starname, since these characters are used to separate the starname from the rest of the filename.

- observations
>   type: string <br>
    default: None <br>
    Full path to the file with the observational data.

- pattern_starting_pulsation
>   type: dictionary of list <br>
    default: None <br>
    Only required when 'provided-pulsation' method is included in 'pattern_methods'. <br>
    This dictionary has two keys, 'period' and 'frequency', holding a list of periods or frequencies, respectively, 
    to start building the pulsation pattern from when using the 'provided-pulsation' method. <br>
    Periods are required when periods or period spacings ('P' or 'dP') are used as observables in `observable_seismic`, 
    and frequencies are required when frequenies ('f') are used as observables in `observable_seismic`. If periods or frequencies are not required, the list can contain `None` instead of values. <br>
    The list should contain one pulsation per part of the split pattern,
    ordered the same as the file with the observations. (List of lenght 1 in case of continuous pattern.)

### Simulated theoretical model grid
- grid_parent_directory
>   type: string <br>
    default: None <br>
    Parent directory of the computed grids.

- grids
>   type: list of strings <br>
    default: None <br>
    Names of the subdirectories of the grids with different physics.

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
    Meridional degree (k value) of the mode ID (k,m) of the g-mode pattern.

- mval
>   type: int <br>
    default: 1 <br>
    Azimuthal order (m value) of the mode ID (k,m) of the g-mode pattern.

- rotation_gyre
>   type: float <br>
    default: 0.0 <br>
    Rotation rate in d^-1 (CYC_PER_DAY of GYRE) that was used in GYRE to calculate the pulsation frequencies.

- gyre_dir
>   type: string <br>
    default: os environment variable $GYRE_DIR <br>
    Path to the GYRE directory.

### Modelling methodology
- pattern_methods
>   type: list of strings <br>
    options: 'provided-pulsation', 'highest-frequency', 'chisq-longest-sequence' <br>
    default: ['provided-pulsation', 'highest-frequency', 'chisq-longest-sequence'] <br>
    List of methods to construct the theoretical frequency pattern (repeats modelling for each method). <br>
    >> 'provided-pulsation' builds the mode pattern starting from the theoretical mode that is closest to the provided pulsation frequency or period. Can be used to provide and start from a trusted mode e.g. the highest-amplitude mode in the observed pattern.
    >
    >> 'highest-frequency' will begin matching modes starting from the theoretical mode that is closest to the highest-frequency mode detected in the observed pattern.
    >
    >>  'chisq-longest-sequence' will match each observed mode period to its best matching theoretical counterpart, and adopt the longest sequence of consecutive modes retrieved in this way. In the case of multiple mode series with the same length, a final pattern selection is made based on the best match between theory and observations. The rest of the theoretical pattern is then build consecutively in radial order starting from this initial 'longest sequence' pattern.

- merit_functions
>   type: list of strings <br>
    options: 'CS', 'MD' <br>
    default: ['CS', 'MD'] <br>
    List of merit functions (repeats modelling for each merit function).
    >> 'CS' uses a reduced chi-squared function.
    >
    >> 'MD' uses Mahalanobis distances, see <a href="https://ui.adsabs.harvard.edu/abs/2018ApJS..237...15A/abstract" target="_blank">Aerts et al. (2018)</a>
    for its application in the context of asteroseismology.

- N_periods
>   type: int <br>
    default: None <br>
    Number of periods (or frequencies) in the observed pattern. Missing frequencies in interrupted patterns are not counted in this number.

- N_pattern_parts
>   type: int <br>
    default: 1 <br>
    In how many parts the observed pattern is split. Defaults to 1 assuming an uninterrupted pattern.

- observable_seismic
>   type: list of strings <br>
    options: 'P', 'dP', 'f' <br>
    default: ['P', 'dP'] <br>
    List of asteroseismic observables to fit in the merit function (repeats modelling for each observable) options are 'P' (period), 'dP' (period-spacing), and 'f' (frequency).

- observable_additional
>   type: list of strings <br>
    default: None <br>
    List of additional observables to use in the merit function.
    Observables other than 'logTeff', 'logL', and 'logg' must correspond to header items included in the MESA profiles.
    Set to None to use only the asteroseismic observables.

- n_sigma_box
>   type: int <br>
    default: 3 <br>
    Ignore models outside of an n-sigma error box on the surface properties (Teff, logL, logg) included in the observational data, set to None to include all models.

- free_parameters
>   type: list <br>
    default: ['Z', 'M', 'logD', 'aov', 'fov', 'Xc'] <br>
    List of varied parameters in the grid, as written in the filenames of the models in the grid (see <a href="{{site.baseurl}}/Walkthrough">walkthrough</a>), that remain free parameters in the modelling.
    In case of a binary system where additional constraints from the companion should be taken into account (see setting `spectro_companion`), this list should start with the first two entries being 'Z', 'M' representing metallicity and mass.

- fixed_parameters
>   type: dict <br>
    default: None <br>
    Dictionary with varied parameters in the grid (dict keys) that are fixed to a certain value (dict value)
    to consider nested grids with fewer free parameters. These parameters should be included in the filenames, just like the free_parameters. Defaults to None to not fix any parameters.
    (E.g. {'aov': 0.0, 'fov': 0.0} to fix both these parameters to zero.)

- evolution_parameter
>   type: string <br>
    default: 'Xc' <br>
    Name of the parameter that is used to track the evolutionary steps of the model.

- evolution_step
>   type: float <br>
    default: -0.01 <br>
    Change in the evolutionary parameter from one step to the next (negative if quantity decreases, e.g. central hydrogen content Xc). Only needed when enforcing binary constraints.

### Modelling binaries and enforcing constraints of the companion star
- constraint_companion
>   type: dictionary <br>
    default: None <br>
    Defaults to None instead of a dictionary to not include any constraints from binarity.
    Dictionary with the following keys holding all binary Information (q = Mass_secondary/Mass_primary): <br> {'q': \<float\>, 'q_err': \<float\>,'Teff': \<float\>, 'Teff_err': \<float\>,
          'logg': \<float\>, 'logg_err': \<float\>, 'logL': \<float\>, 'logL_err':\<float\>, 'primary_pulsates':\<boolean\>}  <br>
    Set one of the observables Teff, logg, logL to None to not use that observable.

- isocloud_grid_directory
>   type: string <br>
    default: None <br>
    The path to the isocloud grid directory. This grid of MESA history files should be constructed in the same way as the regular grid described in <a href="{{site.baseurl}}/Walkthrough">walkthrough</a>, but only needs the history files. This can be the same theoretical grid as the main one that is used, but can also be an additional grid that is sparser and covers a larger range in stellar masses, since the binary companion might not fall within a grid that is more focussed on the pulsating star.

### Other settings
- debugging
>   type: boolean <br>
    default: False <br>
    Set to True to set logger level to debug.

- nr_cpu
>   type:int <br>
    default: None <br>
    Number of worker processes to use in multiprocessing. The default 'None' will cause the pools to use the number returned by os.cpu_count().

- conerplot_axis_labels
>   type: dictionary, keys and values are strings <br>
    default: {'rot': r'$\Omega\_{\mathrm{rot}}$ [d$^{-1}$]', 'M': r'M$\_{\rm ini}$', 'Z': r'Z$\_{\rm ini}$', 'logD':r'log(D$\_{\rm env}$)', 'aov':r'$\alpha\_{\rm CBM}$','fov':r'f$\_{\rm CBM}$','Xc':r'$\rm X_c$'} <br>
    keys are the grid parameters, values are how they should be put in the labels of the cornerplots' axis.

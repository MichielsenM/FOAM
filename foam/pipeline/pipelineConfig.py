import os, sys, logging

class pipelineConfig:
    """
        A python class to keep the configurations settings of the modelling pipeline
    """

    def __init__(self, **kwargs):
        """
        Initialising the instance of the configuration.

        ------- Parameters -------
        -- Info to retrieve observed and simulated data --
        star: string
            Name of the star, used for generating filenames
        observations: string
            Name of (or full path to) the file with the observational data
        periods_or_frequencies_observed: string
            options: 'period' or 'frequency'
            Use the observed periods or frequencies, be consistent in observable_list later.
        highest_amplitude_pulsation: dictionary of lists
            Pulsation with the highest amplitude to build pattern from when using 'highest_amplitude' method.
            List with highest amplitudes per part of the split pattern,
            ordered the same as the file with the observations. (List of lenght 1 in case of continuous pattern.)
        grid_parent_directory: string
            Parent directory of the computed grids
        grids: list of strings
            Names of the subdirectories of the grids with different physics (e.g. different temperature gradients)
        subgrid: string
            String to select a subgrid of the original grid, fixing given parameters

        -- GYRE and mode information --
        kval, mval: int
            Mode ID (k,m) of the g-mode pattern
        rotation_gyre: float
            Rotation rate in d^-1 (CYC_PER_DAY of GYRE)
        gyre_dir: string
            Path to the GYRE directory, uses environment variable GYRE_DIR by default

        -- Modelling methodology --
        pattern_methods: list of strings
            List of methods to construct the theoretical frequency pattern (repeats modelling for each method)
        merit_functions: list of strings
            List of merit functions (repeats modelling for each merit function)
        observable_list: list of list of strings
            Lists of observables to fit (repeats modelling for each list)
            e.g. ['period'] can be expanded to ['period', 'logg'] to include more observables in the likelihood estimation
        observable_aic: list of strings
            calculate AICc for these observables (abbreviated names)
        n_sigma_spectrobox: int
            Ignore models outside of the n-sigma spectroscopic error box, set to None to include all models.

        free_parameters: list
            List of varied parameters in the grid, as written in grid filenames,
            that remain free parameters in the modelling.
        fixed_parameters: dict
            Dictionary with varied parameters in the grid (dict keys) that are fixed to a certain value (dict value)
            to consider nested grids with fewer free parameters. Defaults to None to not fix any parameters.
            (E.g. {'aov': 0.0, 'fov': 0.0} to fix both these mixing parameters to zero.)
        N_periods: int
            Number of periods in the observed pattern
        N_pattern_parts: int
            In how many parts the observed pattern is split. Defaults to 1 assuming an uninterrupted pattern.

        -- For modelling binaries and enfocring constraints of the companion star. --
        spectro_companion: dictionary
            Dictionary with the following keys holding all binary Information (q = Mass_secondary/Mass_primary),
            set a spectro observable to None to not use that observable.
            e.g. {'q': <float>, 'q_err': <float>,'Teff': <float>, 'Teff_err': <float>, 'logg': <float>, 'logg_err': <float>, 'logL': None, 'logL_err':None, 'primary_pulsates':<boolean>}
            Defaults to None instead of a dictionary to not include binary constraints.

        isocloud_grid_directory: string
            path to directory with files containing the isocloud grid summary files
        """
        # Logging settings, other scripts spawn a child logger of this one, copying its settings.
        logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        self.logger = logging.getLogger('logger')
        self.logger.setLevel(logging.DEBUG)


        self.star = kwargs.pop("star", None)
        self.observations = kwargs.pop("observations", None)
        self.periods_or_frequencies_observed = kwargs.pop("periods_or_frequencies_observed", None)
        self.highest_amplitude_pulsation = kwargs.pop("highest_amplitude_pulsation", None)
        self.grid_parent_directory = kwargs.pop("grid_parent_directory", None)
        self.grids = kwargs.pop("grids", None)
        self.subgrid = kwargs.pop("subgrid", "*")

        # Pulsations
        self.kval = kwargs.pop("kval", 0)    # the mode ID of the g-mode pattern
        self.mval = kwargs.pop("mval", 1)
        self.rotation_gyre = kwargs.pop("rotation_gyre", 0.0)
        self.gyre_dir = kwargs.pop("gyre_dir", os.environ.get('GYRE_DIR'))
        if self.gyre_dir is None:
            sys.exit(self.logger.error('Either set GYRE_DIR as environment variable, or specify when initialising the pipelineConfig object'))

        # Modelling methodology
        self.pattern_methods = kwargs.pop("pattern_methods", ['chisq_longest_sequence','highest_amplitude', 'highest_frequency'])
        self.merit_functions = kwargs.pop("merit_functions", ['chi2', 'mahalanobis'])
        self.observable_list = kwargs.pop("observable_list", [['period'], ['period_spacing']])
        self.observable_aic = kwargs.pop("observable_aic", ['P', 'dP'])
        self.n_sigma_spectrobox = kwargs.pop("n_sigma_spectrobox", 3)

        self.free_parameters = kwargs.pop("free_parameters", ['M', 'Z', 'logD', 'aov', 'fov', 'Xc'])
        self.fixed_parameters = kwargs.pop("fixed_parameters", None)
        self.k = len(self.free_parameters) # Number of free parameters
        self.N_periods = kwargs.pop("N_periods", None)
        self.N_pattern_parts = kwargs.pop("N_pattern_parts", 1)
        # Number of observables, calculated from the number of periods
        # Number of period spacings is number of periods minus amount of separated patterns.
        # E.g. uninterrupted pattern: 36 periods, so 35 period spacings
        self.N_dict = {'P' : self.N_periods,'dP': self.N_periods-self.N_pattern_parts}

        # Binarity
        self.spectro_companion = kwargs.pop("spectro_companion", None)
        self.isocloud_grid_directory = kwargs.pop("isocloud_grid_directory", None)

        # Check if all keyword arguments were used
        if len(kwargs) !=0:
            sys.exit(self.logger.error(f'The following keyword arguments were not recognised {[x for x in kwargs.keys()]}'))

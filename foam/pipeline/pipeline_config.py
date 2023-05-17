""" Module for the configuration of the pipeline, contains the class PipelineConfig,
and one instance of this class named 'config'. Through this instance, the configuration
is accessible from all other modules by importing the module and accessing the config attribute of the module."""
import os, sys, logging
from pathlib import Path

class PipelineConfig:
    """
        A python class to keep the configuration settings of the modelling pipeline.
    """

    def __init__(self, **kwargs):
        """
        Initialising the instance of the configuration.

        ------- Parameters -------
        --- Settings about observational data ---
        star: string
            Name of the star, used for generating filenames
        observations: string
            Full path to the file with the observational data
        periods_or_frequencies_observed: string
            options: 'period' or 'frequency'
            Use the observed periods or frequencies, be consistent in observable_list later.
        highest_amplitude_pulsation: dictionary of lists
            Pulsation with the highest amplitude to build pattern from when using 'highest-amplitude' method.
            List with highest amplitudes per part of the split pattern,
            ordered the same as the file with the observations. (List of lenght 1 in case of continuous pattern.)

        --- Simulated theoretical model grid ---
        grid_parent_directory: string
            Parent directory of the computed grids
        grids: list of strings
            Names of the subdirectories of the grids with different physics (e.g. different temperature gradients)
        subgrid: string
            String to select a subgrid of the original grid, fixing given parameters

        --- GYRE and mode information ---
        kval, mval: int
            Mode ID (k,m) of the g-mode pattern
        rotation_gyre: float
            Rotation rate in d^-1 (CYC_PER_DAY of GYRE) that was used in GYRE to calculate the pulsation frequencies
        gyre_dir: string
            Path to the GYRE directory, uses environment variable GYRE_DIR by default

        --- Modelling methodology ---
        pattern_methods: list of strings
            List of methods to construct the theoretical frequency pattern (repeats modelling for each method)
        merit_functions: list of strings
            List of merit functions (repeats modelling for each merit function)
        observable_list: list of list of strings
            Lists of observables to fit (repeats modelling for each list)
            e.g. ['period'] can be expanded to ['period', 'logg'] to include more observables in the likelihood estimation
        observable_aic: list of strings
            calculate AICc for these observables (abbreviated names)
        n_sigma_box: int
            Ignore models outside of the n-sigma error box on the surface properties, set to None to include all models.
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

        --- For modelling binaries and enforcing constraints of the companion star. ---
        constraint_companion: dictionary
            Dictionary with the following keys holding all binary Information (q = Mass_secondary/Mass_primary)
            e.g. {'q': <float>, 'q_err': <float>,'Teff': <float>, 'Teff_err': <float>,
                  'logg': <float>, 'logg_err': <float>, 'logL': <float>, 'logL_err':<float>, 'primary_pulsates':<boolean>}
            set surface observable to None to not use that observable.
            Defaults to None instead of a dictionary to not include any constraints from binarity.

        isocloud_grid_directory: string
            The path to the isocloud grid directory.

        --- Other settings ---
        conerplot_axis_labels: dictionary, its keys and values are strings
            keys are the grid parameters, values are how they should be put in the labels of the cornerplots' axis
        debugging: boolean
            Set to True to set logger level to debug
        nr_cpu: int
            Number of worker processes to use in multiprocessing. The default 'None' will cause the pools to use the number returned by os.cpu_count().
        """
        # Logging settings, other scripts spawn a child logger of this one, copying its settings.
        logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        self.logger = logging.getLogger('logger')
        if kwargs.pop("debugging", False):
            self.logger.setLevel(logging.DEBUG)

        self.nr_cpu = kwargs.pop("nr_cpu", None)

        # Set the main top-level directory
        self.main_directory = os.getcwd()

        # Settings about observational data
        self.star = kwargs.pop("star", None)
        self.observations = kwargs.pop("observations", None)
        self.periods_or_frequencies_observed = kwargs.pop("periods_or_frequencies_observed", 'period')
        self.highest_amplitude_pulsation = kwargs.pop("highest_amplitude_pulsation", None)

        # Simulated theoretical model grid
        self.grid_parent_directory = kwargs.pop("grid_parent_directory", None)
        self.grids = kwargs.pop("grids", None)
        self.subgrid = kwargs.pop("subgrid", "*")

        # Pulsation information
        self.kval = kwargs.pop("kval", 0)    # the mode ID of the g-mode pattern
        self.mval = kwargs.pop("mval", 1)
        self.rotation_gyre = kwargs.pop("rotation_gyre", 0.0)
        self.gyre_dir = kwargs.pop("gyre_dir", os.environ.get('GYRE_DIR'))

        # Modelling methodology
        self.pattern_methods = kwargs.pop("pattern_methods", ['chisq-longest-sequence','highest-amplitude', 'highest-frequency'])
        self.merit_functions = kwargs.pop("merit_functions", ['CS', 'MD'])
        self.observable_list = kwargs.pop("observable_list", [['period'], ['period-spacing']])
        self.observable_aic = kwargs.pop("observable_aic", ['P', 'dP'])
        self.n_sigma_box = kwargs.pop("n_sigma_box", 3)
        self.free_parameters = kwargs.pop("free_parameters", ['Z', 'M', 'logD', 'aov', 'fov', 'Xc'])
        self.fixed_parameters = kwargs.pop("fixed_parameters", None)

        if self.fixed_parameters is None:   # keep track of all relevant parameters, free and fixed, in the grid
            self.grid_parameters = self.free_parameters
        else:
            self.grid_parameters = self.free_parameters+list(self.fixed_parameters.keys())

        self.k = len(self.free_parameters) # Number of free parameters
        self.N_periods = kwargs.pop("N_periods", None)
        self.N_pattern_parts = kwargs.pop("N_pattern_parts", 1)
        # Number of observables, calculated from the number of periods
        # Number of period spacings is number of periods minus amount of separated patterns.
        # E.g. uninterrupted pattern: 36 periods, so 35 period spacings
        self.N_dict = {'P' : self.N_periods,'dP': self.N_periods-self.N_pattern_parts, 'f' : self.N_periods,}

        # Binarity
        self.constraint_companion = kwargs.pop("constraint_companion", None)
        self.isocloud_grid_directory = kwargs.pop("isocloud_grid_directory", None)

        # Plotting
        self.conerplot_axis_labels = kwargs.pop("conerplot_axis_labels", {'rot': r'$\Omega_{\mathrm{rot}}$ [d$^{-1}$]' ,'M': r'M$_{\rm ini}$', 'Z': r'Z$_{\rm ini}$',
                                                'logD':r'log(D$_{\rm env}$)', 'aov':r'$\alpha_{\rm CBM}$','fov':r'f$_{\rm CBM}$','Xc':r'$\rm X_c$'})
        # Check for errors in input arguments
        self._check_init_arguments(kwargs)


    def _check_init_arguments(self, remaining_kwargs):
        """
        Check for errors in the input kwargs. Logs all errors and exits the program if any are found.

        ------- Parameters -------
        remaining_kwargs: list
            kwargs passed to the __init__() that have not been popped by it
        """
        input_error = False

        # Check if all keyword arguments were used
        if len(remaining_kwargs) !=0:
            self.logger.error(f'PipelineConfig: The following keyword arguments were not recognised {[x for x in remaining_kwargs.keys()]}')
            input_error = True

        # Check if a GYRE directory is specified
        if self.gyre_dir is None:
            self.logger.error('PipelineConfig: Either set GYRE_DIR as environment variable, or specify when initialising the PipelineConfig object')
            input_error = True

        # Check star name is provided
        if self.star is None:
            self.logger.error('PipelineConfig: Please provide the name of the modelled star')
            input_error = True

        # Check if file with observations exists
        if self.observations is not None:
            if not Path(self.observations).is_file():
                self.logger.error(f'PipelineConfig: observations file not found: {self.observations} ')
                input_error = True
        else:
            self.logger.error(f'PipelineConfig: observations file not found: {self.observations} ')
            input_error = True

        # Check if parent directory of the grids exists
        if self.grid_parent_directory is not None:
            if not Path(self.grid_parent_directory).is_dir():
                self.logger.error(f'PipelineConfig: Parent directory of theoretical grids not found: {self.grid_parent_directory} ')
                input_error = True

        # Check if directories of grids exists
            for grid in self.grids:
                if not Path(f'{self.grid_parent_directory}/{grid}').is_dir():
                    self.logger.error(f'PipelineConfig: Directory of theoretical grid not found: {self.grid_parent_directory}/{grid} ')
                    input_error = True

        # Check if names of grids are specified
        if self.grids is None:
            self.logger.error(f'PipelineConfig: Name of theoretical grid not specified')
            input_error = True

        # Check that you don't use observed periods whilst looking at the theoretical values as if they are frequencies, and vice versa.
        match_obsAndTheory = False
        for obs_list in self.observable_list:
            for obs in obs_list:
                if (self.periods_or_frequencies_observed) in obs:
                    match_obsAndTheory = True
            if match_obsAndTheory is False:
                self.logger.error(f'The observables that are analysed {self.observable_list} do not all include the observational data that is used: {self.periods_or_frequencies_observed}')
                input_error = True
            match_obsAndTheory = False

        # Check if none of the fixed parameters are in the list of free parameters, and set name for nested grid
        if self.fixed_parameters is not None:
            self.nested_grid_dir = 'Nested_grid_fix'
            for param in self.fixed_parameters.keys():
                self.nested_grid_dir = f'{self.nested_grid_dir}_{param}'
                if param in self.free_parameters:
                    self.logger.error(f'The parameter {param} can not be both fixed and free.')
                    input_error = True

        if input_error:
            sys.exit()
################################################################################

config=None

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
            Path to the file with the observational data
        pattern_starting_pulsation: dictionary of lists
            Dictionary with 'period' and 'frequency' as keys, containing lists of the periods and frequencies, respectively,
            to start building the pattern from when using 'provided-pulsation' method.
            The lists have one pulsation per part of the split pattern,
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
        observable_seismic: list of strings
            List of asteroseismic observables to fit in the merit function (repeats modelling for each observable)
            options are 'P' (period), 'dP' (period-spacing), and 'f' frequency.
        observable_additional: list of strings
            List of additional observables to use in the merit function (e.g. logTeff, logg, logL ...)
            Set to None to just use the asteroseismic observables.
        n_sigma_box: int
            Ignore models outside of the n-sigma error box on the surface properties, set to None to include all models.
        free_parameters: list
            List of varied parameters in the grid, as written in grid filenames, that remain free parameters in the modelling.
            In case of a binary system where additional constraints from the companion should be taken into account, 
            this list should start with the first two entries being 'Z', 'M' representing metallicity and mass.
        fixed_parameters: dict
            Dictionary with varied parameters in the grid (dict keys) that are fixed to a certain value (dict value)
            to consider nested grids with fewer free parameters. Defaults to None to not fix any parameters.
            (E.g. {'aov': 0.0, 'fov': 0.0} to fix both these mixing parameters to zero.)
        evolution_parameter: string
            Name of the parameter that is used to track the evolutionary steps of the model.
        evolution_step: float
            Change in the evolutionary parameter from one step to the next (negative if quantity decreases, e.g. central hydrogen content Xc)       
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
        self.logger.setLevel(logging.INFO)
        if kwargs.pop("debugging", False):
            self.logger.setLevel(logging.DEBUG)

        self.nr_cpu = kwargs.pop("nr_cpu", None)

        # Set the main top-level directory
        self.main_directory = os.getcwd()

        # Settings about observational data
        self.star = kwargs.pop("star", None)
        obs_path = kwargs.pop("observations", None)
        if os.path.isabs(obs_path): # Check if path is absolute path
            self.observations = obs_path
        else:                       # Convert relative path to absolute path
            self.observations = f'{self.main_directory}/{obs_path}'
            
        self.pattern_starting_pulsation = kwargs.pop("pattern_starting_pulsation", None)

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
        self.pattern_methods = kwargs.pop("pattern_methods", ['chisq-longest-sequence','provided-pulsation', 'highest-frequency'])
        self.merit_functions = kwargs.pop("merit_functions", ['CS', 'MD'])
        self.observable_seismic = kwargs.pop("observable_seismic", ['P', 'dP'])
        self.observable_additional = kwargs.pop("observable_additional", None)

        self.n_sigma_box = kwargs.pop("n_sigma_box", 3)
        self.free_parameters = kwargs.pop("free_parameters", ['Z', 'M', 'logD', 'aov', 'fov', 'Xc'])
        self.fixed_parameters = kwargs.pop("fixed_parameters", None)

        if self.fixed_parameters is None:   # keep track of all relevant parameters, free and fixed, in the grid
            self.grid_parameters = self.free_parameters
        else:
            self.grid_parameters = self.free_parameters+list(self.fixed_parameters.keys())

        self.evolution_parameter = kwargs.pop("evolution_parameter", 'Xc')
        self.evolution_step = kwargs.pop("evolution_step", -0.01)
        
        self.k = len(self.free_parameters) # Number of free parameters
        self.N_periods = kwargs.pop("N_periods", None)
        self.N_pattern_parts = kwargs.pop("N_pattern_parts", 1)
        # Number of observables, calculated from the number of periods
        # Number of period spacings is number of periods minus amount of separated patterns.
        # E.g. uninterrupted pattern: 36 periods, so 35 period spacings
        if self.observable_additional is None:
            self.N_dict = {'P' : self.N_periods,'dP': self.N_periods-self.N_pattern_parts, 'f' : self.N_periods}
        else:
            self.N_dict = {'P+extra' : self.N_periods+len(self.observable_additional),'dP+extra': self.N_periods-self.N_pattern_parts+len(self.observable_additional), 'f+extra' : self.N_periods+len(self.observable_additional)}

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
        
        # Check if "Z" and "M" are first free parameters when constraints from binary companion are used.
        if self.constraint_companion is not None:
            if not (self.free_parameters[0] == 'Z' and self.free_parameters[1] == 'M'):
                self.logger.error(f'PipelineConfig: if constraints from a binary companion should be taken into account, "free_parameters" should start with "Z" and "M" as frist entries. '+
                                  f'However the first entries of this list were "{self.free_parameters[0]}" and "{self.free_parameters[1]}".')
                input_error = True

        # Check if the amount of pulsations provided is equal to the amount of parts the pattern is split into.
        if 'provided-pulsation' in self.pattern_methods:
            if ('P' in self.observable_seismic or 'dP' in self.observable_seismic) and not (len(self.pattern_starting_pulsation['period']) == self.N_pattern_parts):
                self.logger.error('To build patterns based on the provided pulsation method, there should be a pulsation provided per part of the (interrupted) pulsation pattern. Incorrect number of periods provided.')
                input_error = True
            if ('f' in self.observable_seismic) and not (len(self.pattern_starting_pulsation['frequency']) == self.N_pattern_parts):
                self.logger.error('To build patterns based on the provided pulsation method, there should be a pulsation provided per part of the (interrupted) pulsation pattern. Incorrect number of frequencies provided.')
                input_error = True

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

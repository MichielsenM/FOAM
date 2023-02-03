"""A configuration file for the pipeline, rename this to 'config.py' and copy to the folder from where you run the pipeline script."""
# Logging settings, other scripts spawn a child logger of this one, copying its settings.
import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('logger')
logger.setLevel(logging.DEBUG)

star = #'KIC7760680'                             # Name of the star
observations = #'data_KIC7760680.tsv'            # File with all the observational data
periods_or_frequencies_observed = 'period'      # # Use the observed periods or frequencies, be consistent in observable_list later
# Pulsation with the highest amplitude to build pattern from in 'highest_amplitude' method. Array with highest amplitude per part of the split pattern
highest_amplitude_pulsation = {'period': [None], 'frequency': [None]} # ordered the same as the file with the observations

gyre_dir = #'/lhome/mathiasm/Software/gyre-6.0.1' # GYRE directory
kval = #0                  # the mode ID of your g-mode pattern
mval = #1
rotation_gyre = # 0.6304

# Parent directory of the computed grid, and names of the directories of the different grids
grid_parent_directory = #'/lhome/mathiasm/MESA_grid_15140'
grids = ['ECP', 'DO']

# Methods to construct the theoretical frequency pattern that are to be used
pattern_methods = ['chisq_longest_sequence','highest_amplitude', 'highest_frequency']

# List of merit functions used
merit_functions = ['chi2', 'mahalanobis']

# Lists of observables to fit
# observable_list = [['period'], ['period_spacing'], ['period', 'logTeff', 'logL', 'logg'], ['frequency']]
observable_list = [['period'], ['period_spacing']]

# calculate AIC for these observables (abbreviated names)
observable_aic = ['P', 'dP']

# String to select a subgrid of the original grid, fixing given parameters
subgrid = '*'
fixed_parameters = None # {'aov': 0.0, 'fov': 0.0}
free_parameters=['M', 'Z', 'logD', 'aov', 'fov', 'Xc']
k = len(free_parameters) # Number of free parameters
N_periods = #36      # Number of periods
N_pattern_parts = 1
# Number of observables, calculated from the number of periods
# Number of period spacings is number of periods minus amount of separated patterns.
# E.g. uninterrupted pattern: 36 periods, so 35 deltaP
N_dict = {'P' : N_periods,'dP': N_periods-N_pattern_parts}

# ignore models outside of the n-sigma spectroscopic error box, set to None to include all models
n_sigma_spectrobox = 3

# For modelling binaries and enfocring constraints of the companion star.
isocloud_grid_directory = None # '/STER/mathiasm/MESA_15140_isochroneGrid'
spectro_companion = None #{'q':0.77, 'q_err':0.09 ,'Teff':14020 , 'Teff_err':280, 'logg':3.55, 'logg_err':0.24, 'logL':None , 'logL_err':None, 'primary_pulsates':True}

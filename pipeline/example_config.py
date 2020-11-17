"""A configuration file for the pipeline, rename this to 'config.py' and copy to the folder from where you run the pipeline script."""
# Logging settings, other scripts spawn a child logger of this one, copying its settings.
import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('logger')
logger.setLevel(logging.DEBUG)

# Name of the star
star = 'KIC7760680'
# File with all the observational data
observations = 'data_KIC7760680.tsv'
periods_or_frequencies_observed = ['period']  # Use the observed periods or frequencies
highest_amplitude_pulsation = {'period': 1.158919, 'frequency': -1}  # pulsation with the highest amplitude to build pattern form in 'highest_amplitude' method

#parent directory of the computed grid, and names of the directories of the different grids
grid_parent_directory = '/lhome/mathiasm/MESA_grid_ECP-DE'
grids = ['ECP', 'DO']

# String to select a subgrid of the original grid, fixing given parameters
subgrid = '*'

# Methods to construct the theoretical frequency pattern
pattern_methods = ['chisq_longest_sequence','highest_amplitude', 'highest_frequency']

# merit functions
merit_functions = ['chi2', 'mahalanobis']

# Observable to fit
# observable_list = [['period'], ['period_spacing'], ['period', 'period_spacing', 'rope_length', 'logTeff', 'logL', 'logg'], ['period', 'period_spacing']]
observable_list = [['period'], ['period_spacing'], ['period', 'period_spacing']]

# calculate AIC for these observables (abbreviated names)
observable_aic = ['P', 'dP', 'P-dP']

# free parameters
k = 6 # M, Z, logD, Xc, aov, fov
free_param=['M', 'Z', 'logD', 'aov', 'fov', 'Xc']
# number of observables, calculated from the number of periods
N_periods = 36
N_dict = {'P' : N_periods,'dP': N_periods-1, 'P-dP': N_periods+N_periods-1}    # 36 periods, so 35 deltaP, so 71 observables

# ignore models outside of the n-sigma spectroscopic error box
n_sigma_spectrobox = 3

################################################################################
# Plotting options
spectroClippedPlots = True

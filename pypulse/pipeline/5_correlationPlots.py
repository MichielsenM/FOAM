"""Make the correlation plots of the grid for the different modelling methodologies."""
from pathlib import Path
import glob
from PyPulse import maximum_likelihood_estimator as mle
from PyPulse import my_python_functions as mypy
import config
################################################################################
if config.spectroClippedPlots:
    files = glob.glob(f'{config.n_sigma_spectrobox}sigmaSpectro_extracted_freqs/*.dat')
else:
    files = glob.glob(f'extracted_freqs/*.dat')

observations = config.observations
for file in files:
    star_name, title = mypy.split_line(Path(file).stem, '_')
    config.logger.info(f'file: {file}')
    mle.plot_correlations(file, observations, fig_title=title, percentile_to_show=0.5, logg_or_logL='logL')

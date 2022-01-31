"""Make the correlation plots of the grid for the different modelling methodologies."""
from pathlib import Path
import glob
from foam import maximum_likelihood_estimator as mle
from foam import support_functions as sf
import config
################################################################################
if config.spectroClippedPlots:
    files = glob.glob(f'{config.n_sigma_spectrobox}sigmaSpectro_extracted_freqs/*.dat')
else:
    files = glob.glob(f'extracted_freqs/*.dat')

observations = config.observations
for file in files:
    # star_name, title = sf.split_line(Path(file).stem, '_')
    title = Path(file).stem
    config.logger.info(f'file: {title}')
    mle.plot_correlations(file, observations, fig_title=title, percentile_to_show=0.5, logg_or_logL='logL')

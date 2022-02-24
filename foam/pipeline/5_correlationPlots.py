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
    Path_file = Path(file)
    title = Path_file.stem
    config.logger.info(f'file: {title}')
    file_ErrorEllips = Path_file.with_stem(f'{Path_file.stem}_error_ellips')
    file_ErrorEllips = str(file_ErrorEllips).replace('extracted_freqs', f'{config.n_sigma_spectrobox}sigmaSpectro_extracted_freqs')
    if not Path(f'figures_correlation/{title}.png').is_file():
        mle.corner_plot(file, file_ErrorEllips, observations, fig_title=title, percentile_to_show=0.5, logg_or_logL='logL')

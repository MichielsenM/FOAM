"""Make the correlation plots of the grid for the different modelling methodologies."""
from pathlib import Path
import glob, multiprocessing
from functools import partial
from foam import plot_tools
from foam.pipeline.pipeline_config import config
################################################################################
if config.n_sigma_spectrobox != None:
    directory_prefix = f'{config.n_sigma_spectrobox}sigmaSpectro_'
else:
    directory_prefix = f''

files = glob.glob(f'extracted_freqs/*[!error_ellips].hdf')

args = []
for file in files:
    Path_file = Path(file)
    title = Path_file.stem
    config.logger.info(f'file: {title}')
    file_ErrorEllips = Path_file.with_stem(f'{Path_file.stem}_2sigma_error_ellipse')
    file_ErrorEllips = str(file_ErrorEllips).replace('extracted_freqs', f'{directory_prefix}extracted_freqs')
    if not Path(file_ErrorEllips).is_file():
        continue
    if not Path(f'{directory_prefix}figures_correlation/{title}.png').is_file():
        args.append((file, file_ErrorEllips, title))

with multiprocessing.Pool() as p:
    func = partial(plot_tools.corner_plot, observations_file=config.observations, fig_outputDir=f'{directory_prefix}figures_correlation/', percentile_to_show=0.5, logg_or_logL='logL', n_sigma_spectrobox=config.n_sigma_spectrobox  )
    p.starmap(func, args)

"""Make the correlation plots of the grid for the different modelling methodologies."""
import matplotlib
from pathlib import Path
import glob, multiprocessing
from functools import partial
from foam import plot_tools
from foam.pipeline.pipeline_config import config
################################################################################
if config.n_sigma_box != None:
    directory_prefix = f'{config.n_sigma_box}sigmaBox_'
else:
    directory_prefix = f''

files = glob.glob(f'meritvalues/*[!error].hdf')

args = []
for file in files:
    Path_file = Path(file)
    title = Path_file.stem
    config.logger.debug(f'file: {title}')
    file_ErrorEllips = Path_file.with_stem(f'{Path_file.stem}_2sigma-error-ellipse')
    file_ErrorEllips = str(file_ErrorEllips).replace('meritvalues', f'{directory_prefix}meritvalues')
    if not Path(file_ErrorEllips).is_file():
        continue
    if not Path(f'{directory_prefix}cornerplots/{title}.png').is_file():
        args.append((file, file_ErrorEllips, title))

matplotlib.use('Agg') # Use this backend for matplotlib to make plots via multiprocessing, otherwise the default gives XIO errors
with multiprocessing.Pool(config.nr_cpu) as p:
    func = partial(plot_tools.corner_plot, observations_file=config.observations, fig_outputDir=f'{directory_prefix}cornerplots/', percentile_to_show=0.5, logg_or_logL='logL', n_sigma_box=config.n_sigma_box, grid_parameters=config.grid_parameters, axis_labels_dict=config.conerplot_axis_labels )
    p.starmap(func, args)

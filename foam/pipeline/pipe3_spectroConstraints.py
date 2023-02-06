"""Ignore all the models that fall outside an n-sigma spectroscopic error box."""
import glob
from pathlib import Path
from foam import maximum_likelihood_estimator as mle
from foam.pipeline.pipelineConfig import config
################################################################################
# Copy of the list of models, and keep only the models that fall within the specified spectroscopic error box
if config.n_sigma_spectrobox != None:
    observations = config.observations
    for grid in config.grids:
        files = glob.glob(f'extracted_freqs/{config.star}_{grid}*.dat')
        for file in files:
            outputFile = f'{config.n_sigma_spectrobox}sigmaSpectro_{file}'
            if not Path(outputFile).is_file():
                mle.spectro_constraint(file, observations, nsigma=config.n_sigma_spectrobox, spectroGrid_file=f'{config.main_directory}/../grid_summary/spectroGrid_{grid}.tsv',
                                    spectro_companion=config.spectro_companion, isocloud_grid_directory=config.isocloud_grid_directory)

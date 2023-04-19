"""Ignore all the models that fall outside an n-sigma spectroscopic error box."""
import glob
from pathlib import Path
import pandas as pd
import multiprocessing
from functools import partial
from foam import additional_constraints as ac
from foam import model_grid as mg
from foam.pipeline.pipeline_config import config

################################################################################
# Copy of the list of models, and keep only the models that fall within the specified spectroscopic error box
if config.n_sigma_spectrobox != None:
    observations = config.observations

    params = list(config.free_parameters) # To make a copy and not remove Xc from the config
    params.remove('Xc')
    summary = mg.GridSummary(params)

    if not Path('isocloud_grid.h5').is_file():
        summary.create_summary_file(config.isocloud_grid_directory, columns=['star_age','log_L','log_Teff','log_g'], magnitudes=False, output_name='isocloud_grid.h5', file_ending='hist', files_directory_name='history')
    else:
        summary.read_summary_file('isocloud_grid.h5')

    isocloud_summary_dict = {}
    for Z in summary.Z_array:
        isocloud_summary_dict.update({Z:{}})
    for Z in summary.Z_array:
        for M in summary.M_array:
            df_MZ = pd.DataFrame()
            for logD in summary.logD_array:
                for aov in summary.aov_array:
                    for fov in summary.fov_array:

                        df = pd.DataFrame( data = summary.grid_data[f'{Z:.3f}'][f'{M:.2f}'][f'{logD:.2f}'][f'{aov:.3f}'][f'{fov:.3f}'])
                        df_MZ = pd.concat([df_MZ, df], ignore_index=True)
            isocloud_summary_dict[Z].update({M : df_MZ })


    all_files_kept = []
    for grid in config.grids:
        files = glob.glob(f'extracted_freqs/{config.star}_{grid}*.hdf')
        files_kept = list(files)
        for file in files:
            output_file = f'{config.n_sigma_spectrobox}sigmaSpectro_{file}'
            if Path(output_file).is_file():
                files_kept.remove(file)
                config.logger.warning(f'file already existed: {output_file}')
        all_files_kept.extend(files_kept)

    with multiprocessing.Pool(4) as p: # For some reason 4 is fastest, find out why more becomes slower
        func = partial( ac.spectro_constraint,  observations_file=observations, nsigma=config.n_sigma_spectrobox, spectroGrid_file=f'{config.main_directory}/../grid_summary/spectroGrid_{grid}.hdf',
                            spectro_companion=config.spectro_companion, isocloud_grid_summary=isocloud_summary_dict)
        p.map(func, all_files_kept)

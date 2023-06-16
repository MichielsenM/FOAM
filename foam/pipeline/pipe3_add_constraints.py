"""Ignore all the models that fall outside an n-sigma error box."""
import glob
from pathlib import Path
import pandas as pd
import multiprocessing
from functools import partial
from foam import additional_constraints as ac
from foam import model_grid as mg
from foam.pipeline.pipeline_config import config

################################################################################
# Copy of the list of models, and keep only the models that fall within the specified error box
if config.n_sigma_box != None:
    observations = config.observations

    isocloud_summary_dict = None
    if config.constraint_companion is not None:
        if not Path(f'{config.main_directory}/isocloud_grid.h5').is_file():
            params = list(config.free_parameters) # To make a copy and not remove Xc from the config
            params.extend(config.fixed_parameters)
            params.remove('Xc')
            summary = mg.GridSummary(params)
            summary.create_summary_file(config.isocloud_grid_directory, columns=['star_age','log_L','log_Teff','log_g'], magnitudes=False, output_name=f'{config.main_directory}/isocloud_grid.h5', file_ending='hist', files_directory_name='history')
        else:
            summary = mg.GridSummary(None)
            summary.read_summary_file(f'{config.main_directory}/isocloud_grid.h5')

        # Create dictionary with all the 'star_age','log_L','log_Teff','log_g' values of the whole isocloud per mass-metallicity combination
        isocloud_summary_dict = {}
        for Z in summary.Z_array:
            isocloud_summary_dict.update({Z:{}})
        for Z in summary.Z_array:
            for M in summary.M_array:
                df_MZ = pd.DataFrame()
                for logD in summary.logD_array:
                    for aov in summary.aov_array:
                        for fov in summary.fov_array:
                            df = pd.DataFrame( data = summary.grid_data[f'{Z}'][f'{M}'][f'{logD}'][f'{aov}'][f'{fov}'])
                            df_MZ = pd.concat([df_MZ, df], ignore_index=True)
                isocloud_summary_dict[Z].update({M : df_MZ })


    files_to_analyse = []
    for grid in config.grids:
        files = glob.glob(f'meritvalues/{config.star}_{grid}*.hdf')
        files_kept = list(files)
        for file in files:
            output_file = f'{config.n_sigma_box}sigmaBox_{file}'
            if Path(output_file).is_file():
                files_kept.remove(file)
                config.logger.warning(f'file already existed: {output_file}')
        files_to_analyse.extend(files_kept)

    nr_cpu = 4
    if config.nr_cpu is not None:
        nr_cpu = min(config.nr_cpu, 4)
    with multiprocessing.Pool( nr_cpu ) as p: # For some reason 4 processes is faster than more, find out why more becomes slower
        func = partial( ac.surface_constraint,  observations_file=observations, nsigma=config.n_sigma_box, surfaceGrid_file=f'{config.main_directory}/../grid_summary/surfaceGrid_{grid}.hdf',
                            constraint_companion=config.constraint_companion, isocloud_grid_summary=isocloud_summary_dict, free_parameters = config.free_parameters, evolution_parameter=config.evolution_parameter, evolution_step = config.evolution_step)
        p.map(func, files_to_analyse)

"""Ignore all the models that fall outside an n-sigma spectroscopic error box."""
import glob
from pathlib import Path
import pandas as pd
import multiprocessing
from functools import partial
from foam import additional_constraints as ac
from foam import modelGrid as mg
from foam.pipeline.pipelineConfig import config

################################################################################
# Copy of the list of models, and keep only the models that fall within the specified spectroscopic error box
if config.n_sigma_spectrobox != None:
    observations = config.observations

    params = list(config.free_parameters) # To make a copy and not remove Xc from the config
    params.remove('Xc')
    summary = mg.gridSummary(params)

    import time
    t1=time.perf_counter()
    if not Path('isocloud_grid.h5').is_file():
        summary.create_summary_file(config.isocloud_grid_directory, columns=['star_age','log_L','log_Teff','log_g'], magnitudes=False, output_name='isocloud_grid.h5', file_ending='hist', files_directory_name='history')
        t2=time.perf_counter()
        print(f'Creating isocloud took: {t2-t1}')
    else:
        summary.read_summary_file('isocloud_grid.h5')

        t2=time.perf_counter()
        print(f'Reading isocloud took: {t2-t1}')

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

    t3=time.perf_counter()
    print(f'Making dicts took: {t3-t2}')

    for grid in config.grids:
        files = glob.glob(f'extracted_freqs/{config.star}_{grid}*.hdf')
        files_kept = list(files)
        for file in files:
            output_file = f'{config.n_sigma_spectrobox}sigmaSpectro_{file}'
            if Path(output_file).is_file():
                files_kept.remove(file)
                config.logger.warning(f'file already existed: {output_file}')

        with multiprocessing.Pool(1) as p:
            func = partial( ac.spectro_constraint,  observations_file=observations, nsigma=config.n_sigma_spectrobox, spectroGrid_file=f'{config.main_directory}/../grid_summary/spectroGrid_{grid}.hdf',
                                spectro_companion=config.spectro_companion, isocloud_grid_summary=isocloud_summary_dict)
            # for result in p.map(func, files_kept):
                # item=result
            # p.map_async(func, files_kept)
            p.map(func, files_kept)
            p.close()
            p.join()


        # for file in files:
        #     output_file = f'{config.n_sigma_spectrobox}sigmaSpectro_{file}'
        #     if not (Path(output_file).is_file()):
        #         ac.spectro_constraint(file, observations_file=observations, nsigma=config.n_sigma_spectrobox, spectroGrid_file=f'{config.main_directory}/../grid_summary/spectroGrid_{grid}.hdf',
        #                             spectro_companion=config.spectro_companion, isocloud_grid_summary=isocloud_summary_dict)
        #     else:
        #         config.logger.warning(f'file already existed: {output_file}')

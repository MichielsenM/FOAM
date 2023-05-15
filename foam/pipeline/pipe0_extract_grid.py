""" From the grid location, extract all spectroscopic info from MESA profiles and all frequencies from GYRE models."""
from pathlib import Path
from foam import functions_for_gyre as ffg
from foam import functions_for_mesa as ffm
from foam.pipeline.pipeline_config import config
################################################################################

for grid in config.grids:
    output_file = f'../grid_summary/spectroGrid_{grid}.hdf'
    if not Path(output_file).is_file():
        ffm.grid_extract_spectroscopy(f'{config.grid_parent_directory}/{grid}/MESA_out/Zini*/profiles/*{config.subgrid}*prof', output_file=output_file)
    else:
        config.logger.warning(f'file already existed: {output_file}')

    output_file = f'../grid_summary/pulsationGrid_{grid}_rot{config.rotation_gyre}_k{config.kval}m{config.mval}.hdf'
    if not Path(output_file).is_file():
        ffg.extract_frequency_grid(f'{config.grid_parent_directory}/{grid}/GYRE_out/rot{config.rotation_gyre}_k{config.kval}m{config.mval}/Zini*/*{config.subgrid}*.HDF', output_file=output_file)
    else:
        config.logger.warning(f'file already existed: {output_file}')

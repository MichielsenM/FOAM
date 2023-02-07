""" From the grid location, extract all spectroscopic info from MESA profiles and all frequencies from GYRE models."""
from foam import functions_for_gyre as ffg
from foam import functions_for_mesa as ffm
from foam.pipeline.pipelineConfig import config
################################################################################

if config.mval>0:
    rotation_spin_direction = 'prograde'
elif config.mval==0:
    rotation_spin_direction = 'zonal'
elif config.mval<0:
    rotation_spin_direction = 'retrograde'

for grid in config.grids:
    output_file = f'../grid_summary/spectroGrid_{grid}.hdf'
    ffm.grid_extract_spectroscopy(f'{config.grid_parent_directory}/{grid}/MESA_out/Zini*/profiles/*{config.subgrid}*prof', output_file=output_file)

    output_file = f'../grid_summary/pulsationGrid_{grid}_{config.rotation_gyre}_{rotation_spin_direction}.hdf'
    ffg.extract_frequency_grid(f'{config.grid_parent_directory}/{grid}/GYRE_out/rot{config.rotation_gyre}_{rotation_spin_direction}/Zini*/*{config.subgrid}*.HDF', output_file=output_file)

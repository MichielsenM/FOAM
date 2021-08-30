""" From the grid location, extract all spectroscopic info from MESA profiles and all frequencies from GYRE models."""
from foam import functions_for_gyre as ffg
from foam import functions_for_mesa as ffm
import config # imports the config file relative to the location of the main script
################################################################################
for grid in config.grids:
    output_file = f'grid_summary/spectroGrid_{grid}.tsv'
    ffm.grid_extract_spectroscopy(f'{config.grid_parent_directory}/{grid}/MESA_out/Zini*/profiles/*{config.subgrid}*prof', output_file=output_file)

    output_file = f'grid_summary/pulsationGrid_{grid}.tsv'
    ffg.extract_frequency_grid(f'{config.grid_parent_directory}/{grid}/GYRE_out/rot*/Zini*/*{config.subgrid}*.HDF', output_file=output_file)

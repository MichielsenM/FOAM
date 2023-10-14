"""Construct the theoretical pulsation patterns and merge with surface info into one file."""
from pathlib import Path
from foam import support_functions as sf
from foam import gmode_rotation_scaling as grs
from foam import build_optimised_pattern as bop
from foam.pipeline.pipeline_config import config
################################################################################

# The required asymptotic class object
asymp_obj = grs.Asymptotic(gyre_dir=config.gyre_dir, kval=config.kval, mval=config.mval)

# Construct the pulsation patterns according to the different methods for the extracted theoretical grids
for grid in config.grids:
    surface= f'../grid_summary/surfaceGrid_{grid}.hdf'
    for method in config.pattern_methods:   # Methods to construct theoretical pulsation patterns
        observed_quantities = []            # Whether periods or frequenies are used
        if 'P' or 'dP' in config.observable_seismic:
            observed_quantities.append('period')
        if 'f' in config.observable_seismic:
            observed_quantities.append('frequency')

        for observed_quantity in observed_quantities:
            puls_file = f'extracted_freqs/{observed_quantity}_{config.star}_{grid}_{method}.hdf'
            if not Path(puls_file).is_file():
                bop.construct_theoretical_puls_pattern(f'../grid_summary/pulsationGrid_{grid}_rot{config.rotation_gyre}_k{config.kval}m{config.mval}.hdf', config.observations,
                                                    method, pattern_starting_pulsation=config.pattern_starting_pulsation[observed_quantity], which_observable=observed_quantity, output_file=puls_file,
                                                    asymptotic_object=asymp_obj, estimated_rotation=config.rotation_gyre, grid_parameters=config.grid_parameters, nr_cpu=config.nr_cpu)
            else:
                config.logger.warning(f'file already existed: {puls_file}')

            # Merge surface and pulsation info into one file
            output_name = f'{Path(puls_file).parent}/surface+{Path(puls_file).name}'
            if not Path(output_name).is_file():
                sf.add_surface_to_puls_grid(puls_file, surface, output_name, grid_parameters=config.grid_parameters)
            else:
                config.logger.warning(f'file already existed: {output_name}')

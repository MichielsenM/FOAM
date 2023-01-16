"""Construct the theoretical pulsation patterns and merge with spectroscopic info into one file."""
from pathlib import Path
from foam import functions_for_mesa as ffm
from foam import gmode_rotation_scaling as grs
from foam import build_optimised_pattern as bop
from foam.pipeline.pipelineConfig import config
################################################################################

# The required asymptotic class object
asymp_obj = grs.asymptotic(gyre_dir=config.gyre_dir, kval=config.kval, mval=config.mval)

if config.mval>0:
    rotation_spin_direction = 'prograde'
elif config.mval==0:
    rotation_spin_direction = 'zonal'
elif config.mval<0:
    rotation_spin_direction = 'retrograde'

# Construct the pulsation patterns according to the different methods for the extracted theoretical grids
observable = config.periods_or_frequencies_observed
for grid in config.grids:
    spectro= f'../grid_summary/spectroGrid_{grid}.tsv'
    for method in config.pattern_methods:   # Methods to construct theoretical pulsation patterns
        puls_file = f'extracted_freqs/{observable}_{config.star}_{grid}_{method}.tsv'
        bop.construct_theoretical_freq_pattern(f'../grid_summary/pulsationGrid_{grid}_{config.rotation_gyre}_{rotation_spin_direction}.tsv', config.observations, method, highest_amplitude_pulsation=config.highest_amplitude_pulsation[observable], which_observable=observable, output_file=puls_file, asymptotic_object=asymp_obj, estimated_rotation=config.rotation_gyre)

        # Merge spectro and pulsation info into one file
        output_name = f'{Path(puls_file).parent}/spectro+{Path(puls_file).name}'
        ffm.add_spectro_to_puls_grid(puls_file, spectro, output_name)
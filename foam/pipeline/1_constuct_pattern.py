"""Construct the theoretical pulsation patterns and merge with spectroscopic info into one file."""
from pathlib import Path
from foam import functions_for_gyre as ffg
from foam import functions_for_mesa as ffm
import config # imports the config file relative to the location of the main script
################################################################################
# Construct the pulsation patterns according to the different methods for the extracted theoretical grids
observable = config.periods_or_frequencies_observed
for grid in config.grids:
    spectro= f'grid_summary/spectroGrid_{grid}.tsv'
    for method in config.pattern_methods:   # Methods to construct theoretical pulsation patterns
        puls_file = f'extracted_freqs/{observable}_{config.star}_{grid}_{method}.tsv'
        ffg.construct_theoretical_freq_pattern(f'grid_summary/pulsationGrid_{grid}_{config.rotation}.tsv', config.observations, method, highest_amplitude_pulsation=config.highest_amplitude_pulsation[observable], which_observable=observable, output_file=puls_file)

        # Merge spectro and pulsation info into one file
        output_name = f'{Path(puls_file).parent}/spectro+{Path(puls_file).name}'
        ffm.add_spectro_to_puls_grid(puls_file, spectro, output_name)

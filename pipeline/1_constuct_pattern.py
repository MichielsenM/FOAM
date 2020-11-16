"""Construct the theoretical pulsation patterns and merge with spectroscopic info into one file."""
import glob
from pathlib import Path
from PyPulse import functions_for_gyre as ffg
from PyPulse import functions_for_mesa as ffm
import config # imports the config file relative to the location of the main script
################################################################################
# Construct the pulsation patterns according to the different methods for the extracted theoretical grids
for observable in config.periods_or_frequencies_observed:
    for grid in config.grids:
        for method in config.pattern_methods:   # Methods to construct theoretical pulsation patterns
            out_file = f'extracted_freqs/{observable}_{config.star}_{grid}_{method}.tsv'
            ffg.construct_theoretical_freq_pattern(f'grid_summary/pulsationGrid_{grid}.tsv', config.observations, method, highest_amplitude_pulsation=config.highest_amplitude_pulsation[observable], which_observable=observable, output_file=out_file)
################################################################################
# Merge spectro and pulsation info into one file
for grid in config.grids:
    spectro= f'grid_summary/spectroGrid_{grid}.tsv'
    freq_files = glob.glob(f'extracted_freqs/*_{grid}_*.tsv')
    for freqs in freq_files:
        output_name = f'{Path(freqs).parent}/spectro+{Path(freqs).name}'
        ffm.add_spectro_to_puls_grid(freqs, spectro, output_name)

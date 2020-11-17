"""Calculate the likelihood of all the theoretial patterns according to the specified merit functions."""
from PyPulse import maximum_likelihood_estimator as mle
from functools import partial
import multiprocessing
from pathlib import Path
import config # imports the config file relative to the location of the main script
###############################################################################
file_Path = Path(f'V_matrix/determinant_conditionNr.tsv')
if file_Path.is_file(): file_Path.unlink()  #remove file if it exists to avoid duplicate entries on successive runs
p = multiprocessing.Pool()	# Multiprocessing pool, NR CORES = nr items in obeservables list
for grid in config.grids:
    for method in config.pattern_methods:
        for merit_function in config.merit_functions:
            Theo_path = f'extracted_freqs/spectro+period_{config.star}_{grid}_{method}.tsv'

            func = partial(mle.calculate_likelihood, config.observations, Theo_path, merit_function = merit_function)
            for result in p.imap(func, config.observable_list):
                item=result

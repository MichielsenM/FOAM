"""Calculate the likelihood of all the theoretical patterns according to the specified merit functions."""
from foam import maximum_likelihood_estimator as mle
from functools import partial
import multiprocessing
from pathlib import Path
from foam.pipeline.pipelineConfig import config
###############################################################################
file_Path = Path(f'V_matrix/{config.star}_determinant_conditionNr.tsv')
if file_Path.is_file(): file_Path.unlink()  #remove file if it exists to avoid duplicate entries on successive runs
p = multiprocessing.Pool()	# Multiprocessing pool, NR CORES = nr items in obeservables list
for grid in config.grids:
    for method in config.pattern_methods:
        for merit_function in config.merit_functions:
            Theo_path = f'{config.main_directory}/extracted_freqs/spectro+{config.periods_or_frequencies_observed}_{config.star}_{grid}_{method}.hdf'

            func = partial(mle.calculate_likelihood, config.observations, Theo_path, merit_function = merit_function, star_name=config.star, fixed_params=config.fixed_parameters)
            for result in p.imap(func, config.observable_list):
                item=result
p.close()

"""Calculate the likelihood of all the theoretical patterns according to the specified merit functions."""
from foam import maximum_likelihood_estimator as mle
from functools import partial
import multiprocessing
from pathlib import Path
from foam.pipeline.pipeline_config import config
###############################################################################
file_Path = Path(f'V_matrix/{config.star}_determinant_conditionNr.tsv')
if file_Path.is_file(): file_Path.unlink()  #remove file if it exists to avoid duplicate entries on successive runs
args = []
for grid in config.grids:
    for method in config.pattern_methods:
        Theo_path = f'{config.main_directory}/extracted_freqs/surface+{config.periods_or_frequencies_observed}_{config.star}_{grid}_{method}.hdf'
        for merit_function in config.merit_functions:
            for obs in config.observable_list:
                args.append(( Theo_path, obs, merit_function ))

with multiprocessing.Pool(config.nr_cpu) as p:
    func = partial(mle.calculate_likelihood, Obs_path=config.observations, star_name=config.star, fixed_params=config.fixed_parameters, grid_parameters=config.grid_parameters)
    p.starmap(func, args)

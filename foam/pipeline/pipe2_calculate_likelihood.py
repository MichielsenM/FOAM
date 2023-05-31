"""Calculate the likelihood of all the theoretical patterns according to the specified merit functions."""
from foam import maximum_likelihood_estimator as mle
from functools import partial
import multiprocessing, os
from pathlib import Path
from foam.pipeline.pipeline_config import config
###############################################################################
DataOutDir = Path(f'{os.getcwd()}/meritvalues')
Path(DataOutDir).mkdir(parents=True, exist_ok=True)

file_Path = Path(f'V_matrix/{config.star}_determinant_conditionNr.tsv')
if file_Path.is_file(): file_Path.unlink()  #remove file if it exists to avoid duplicate entries on successive runs
args = []
observables = []
for grid in config.grids:
    for method in config.pattern_methods:
        for merit_function in config.merit_functions:
            for obs in config.observable_seismic:
                if obs=='P' or obs=='dP':
                    observed_quantity = 'period'
                elif obs == 'f':
                    observed_quantity = 'frequency'

                Theo_path = f'{config.main_directory}/extracted_freqs/surface+{observed_quantity}_{config.star}_{grid}_{method}.hdf'
                observables = [obs]
                if config.observable_additional is not None:
                    observables += config.observable_additional
                args.append(( Theo_path, observables, merit_function ))

with multiprocessing.Pool(config.nr_cpu) as p:
    func = partial(mle.calculate_likelihood, Obs_path=config.observations, star_name=config.star, fixed_params=config.fixed_parameters, grid_parameters=config.grid_parameters)
    p.starmap(func, args)

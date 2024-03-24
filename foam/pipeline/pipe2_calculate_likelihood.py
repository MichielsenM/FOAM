"""Calculate the likelihood of all the theoretical patterns according to the specified merit functions."""

import multiprocessing
import os
from functools import partial
from pathlib import Path

from foam import maximum_likelihood_estimator as mle
from foam.pipeline.pipeline_config import config

###############################################################################
data_out_dir = Path(f"{os.getcwd()}/meritvalues")
Path(data_out_dir).mkdir(parents=True, exist_ok=True)

file_path = Path(f"V_matrix/{config.star}_determinant_conditionNr.tsv")
# remove file if it exists to avoid duplicate entries on successive runs
if file_path.is_file():
    file_path.unlink()
args = []
observables = []
for grid in config.grids:
    for method in config.pattern_methods:
        for merit_function in config.merit_functions:
            for obs in config.observable_seismic:
                if obs == "P" or obs == "dP":
                    observed_quantity = "period"
                elif obs == "f":
                    observed_quantity = "frequency"

                theory_path = f"{config.main_directory}/extracted_freqs/surface+{observed_quantity}_{config.star}_{grid}_{method}.hdf"
                observables = [obs]
                if config.observable_additional is not None:
                    observables += config.observable_additional
                args.append((theory_path, observables, merit_function))

with multiprocessing.Pool(config.nr_cpu) as p:
    func = partial(
        mle.calculate_likelihood,
        obs_path=config.observations,
        star_name=config.star,
        fixed_params=config.fixed_parameters,
        grid_parameters=config.grid_parameters,
    )
    p.starmap(func, args)

"""Write the best models of the grid as a LaTeX table."""
import pandas as pd
from foam.pipeline.pipeline_config import config
################################################################################
if config.n_sigma_box != None:
    directory_prefix = f'{config.n_sigma_box}sigmaBox_'
else:
    directory_prefix = f''

if config.observable_additional is not None:
    extra_obs = '+extra'
else:
    extra_obs = ''
################################################################################

for merit in config.merit_functions:
    # Get the pre-calculated AICc values from another file
    df_AICc = pd.read_table(f'{directory_prefix}output_tables/{config.star}_AICc-values_{merit}.tsv', delim_whitespace=True, header=0)
    with open(f'{directory_prefix}output_tables/{config.star}_best-model-table_{merit}.txt', 'w') as outfile:
        params = ''
        for p in config.grid_parameters:
            params += f' {p}' 
        outfile.write(f'Grid Observables Pattern_construction{params} Omega_rot {merit} AICc_{merit}\n')

    best_model_dict = {}
    # Make a dictionary of the best models for each grid and observable combo, according to each merit function
    for pattern in config.pattern_methods:
        for grid in config.grids:
            for obs in config.observable_seismic:
                obs+=extra_obs
                MLE_values_file = f'{directory_prefix}meritvalues/{config.star}_{grid}_{pattern}_{merit}_{obs}_2sigma-error-ellipse.hdf'
                df = pd.read_hdf(MLE_values_file)
                best_model = df.loc[df['meritValue'].idxmin()]

                best_model_dict.update({f'{grid} {merit} {obs} {pattern}': best_model})

                line = f'{grid} {obs} {pattern}'
                for p in config.grid_parameters:
                    line += f' {best_model[p]}'
                line += f' {round(best_model["rot"], 4)}'
                line += f' {round(best_model["meritValue"], 2)}'
                
                name = f'{config.star}_{grid}_{pattern}_{merit}_{obs}'
                line += f' {round(df_AICc.loc[df_AICc.method == name, "AICc"].values[0], 2)}'

                with open(f'{directory_prefix}output_tables/{config.star}_best-model-table_{merit}.txt', 'a') as outfile:
                    outfile.write(f'{line}\n')
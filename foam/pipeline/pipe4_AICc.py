""" Calculate the AICc (Akaike information criterion, corrected for small sample size) and write to a tsv. """
import glob, sys
import pandas as pd
import numpy as np
from pathlib import Path
from foam import support_functions as sf
from foam.pipeline.pipeline_config import config
################################################################################
k = config.k            # number of free paramters in the grid
if config.n_sigma_box != None:
    directory_prefix = f'{config.n_sigma_box}sigmaBox_'
else:
    directory_prefix = f''

if config.observable_additional is not None:
    extra_obs = '+extra'
else:
    extra_obs = ''

output_folder = f'{directory_prefix}output_tables'
Path(output_folder).mkdir(parents=True, exist_ok=True)
################################################################################
# Get the condition numbers file to use its listed values of ln(det(V)) with V the variance-covariance matrix.
for merit in config.merit_functions:
    if merit == 'CS':
        df_AICc = pd.DataFrame( data = [], columns=['method', 'AICc'] )
    elif merit == 'MD':
        df_AICc = pd.read_table(f'V_matrix/{config.star}_determinant_conditionNr.tsv', delim_whitespace=True, header=0)

    for grid in config.grids:
        for method in config.pattern_methods:
            for obs in config.observable_seismic:
                obs+=extra_obs
                df = pd.read_hdf(f'{directory_prefix}meritvalues/{config.star}_{grid}_{method}_{merit}_{obs}.hdf')
                df = df.sort_values('meritValue', ascending=True)

                # Calculate the AICc
                N = config.N_dict[obs] # number of observables
                if merit == 'CS':
                    AICc = (df["meritValue"].iloc[0]/(N-k)) + (2*k*N)/(N-k-1)
                    df_AICc.loc[len(df_AICc)] = [ f'{config.star}_{grid}_{method}_CS_{obs}' , AICc]

                elif merit == 'MD':
                    lndetV =  df_AICc.loc[df_AICc['method'] == f'{config.star}_{grid}_{method}_MD_{obs}', 'ln(det(V))']
                    AICc = df["meritValue"].iloc[0] + k*np.log(2*np.pi) + lndetV + (2*k*N)/(N-k-1)
                    df_AICc.loc[df_AICc.method == f'{config.star}_{grid}_{method}_MD_{obs}', 'AICc'] = AICc

    df_AICc.to_csv(f'{output_folder}/{config.star}_AICc-values_{merit}.tsv', sep='\t',index=False)

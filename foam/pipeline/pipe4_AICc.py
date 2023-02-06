""" Calculate the AICc (Akaike information criterion, corrected for small sample size) and write to a tsv. """
import glob, sys
import pandas as pd
import numpy as np
from pathlib import Path
from foam import support_functions as sf
from foam.pipeline.pipelineConfig import config
################################################################################
k = config.k            # number of free paramters in the grid
merit_abbrev = {'chi2': 'CS', 'mahalanobis': 'MD'}
if config.n_sigma_spectrobox != None:
    directory_prefix = f'{config.n_sigma_spectrobox}sigmaSpectro_'
else:
    directory_prefix = f''
################################################################################
# Get the condition numbers file to use its listed values of ln(det(V)) with V the variance-covariance matrix.
condition_nr_file = f'V_matrix/{config.star}_determinant_conditionNr.tsv'
try:
    df_AICc_MD = pd.read_table(condition_nr_file, delim_whitespace=True, header=0)
except:
    sys.exit(config.logger.error(f'File does not exist: {condition_nr_file} \n rerun the "calculate_likelihood" script using Mahalanobis Distances, \
                            \n or make a file with the same "Method" column (including column entries) in case you just want to use chi-squared.'))

df_AICc_Chi2 = df_AICc_MD[['method']].copy()
for i in range(0, df_AICc_Chi2['method'].size):
    df_AICc_Chi2.iloc[i]['method'] = df_AICc_Chi2.iloc[i]['method'].replace('MD', 'CS')

for merit in config.merit_functions:
    merit = merit_abbrev[merit]
    for obs in config.observable_aic:
        files = glob.glob(f'{directory_prefix}extracted_freqs/*{merit}_{obs}.dat')
        for file in sorted(files):
            Path_file = Path(file)
            star_name, analysis = sf.split_line(Path_file.stem, '_')
            df = pd.read_csv(file, delim_whitespace=True, header=0)
            df = df.sort_values('meritValue', ascending=True)

            # Calculate the AICc
            N = config.N_dict[obs] # number of observables
            if merit == 'CS':
                AICc = (df["meritValue"].iloc[0]/(N-k)) + (2*k*N)/(N-k-1)
                df_AICc_Chi2.loc[df_AICc_Chi2.method == f'{config.star}_{analysis}', 'AICc'] = AICc

            elif merit == 'MD':
                lndetV =  df_AICc_MD.loc[df_AICc_MD['method'] == f'{config.star}_{analysis}', 'ln(det(V))']
                AICc = df["meritValue"].iloc[0] + k*np.log(2*np.pi) + lndetV + (2*k*N)/(N-k-1)
                df_AICc_MD.loc[df_AICc_MD.method == f'{config.star}_{analysis}', 'AICc'] = AICc

output_folder = f'{directory_prefix}output_tables'
Path(output_folder).mkdir(parents=True, exist_ok=True)
df_AICc_MD.to_csv(f'{output_folder}/{config.star}_AICc_values_MD.tsv', sep='\t',index=False)
df_AICc_Chi2.to_csv(f'{output_folder}/{config.star}_AICc_values_Chi2.tsv', sep='\t',index=False)

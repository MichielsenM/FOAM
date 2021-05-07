"""Ignore all the models that fall outside an n-sigma spectroscopic error box.
Calculate the AICc (Akaike information criterion, corrected for small sample size) and write to a tsv.
"""
import glob, sys
import pandas as pd
import numpy as np
from pathlib import Path
from PyPulse import my_python_functions as mypy
from PyPulse import maximum_likelihood_estimator as mle
import config # imports the config file relative to the location of the main script
################################################################################
# Copy of the list of models, and keep only the models that fall within the specified spectroscopic error box
files = glob.glob(f'extracted_freqs/*.dat')
observations = config.observations
for file in files:
    mle.spectro_cutoff(file, observations, nsigma=config.n_sigma_spectrobox)
################################################################################
k = config.k            # number of free paramters in the grid
merit_abbrev = {'chi2': 'CS', 'mahalanobis': 'MD'}
################################################################################
# Get the condition numbers file to use its listed values of ln(det(V)) with V the variance-covariance matrix.
condition_nr_file = 'V_matrix/determinant_conditionNr.tsv'
try:
    df_AICc_MD = pd.read_table(condition_nr_file, delim_whitespace=True, header=0)
except:
    sys.exit(config.logger.error(f'File does not exist: {condition_nr_file} \n rerun the "calculate_likelihood" script using Mahalanobis Distances, \
                            \n or make a file with the same "Method" column (including column entries) in case you just want to use chi-squared.'))

df_AICc_Chi2 = df_AICc_MD[['method']].copy()
for i in range(0, df_AICc_Chi2.size):
    df_AICc_Chi2.iloc[i]['method'] = df_AICc_Chi2.iloc[i]['method'].replace('MD', 'CS')

for merit in config.merit_functions:
    merit = merit_abbrev[merit]
    for obs in config.observable_aic:
        files = glob.glob(f'{config.n_sigma_spectrobox}sigmaSpectro_extracted_freqs/*{merit}_{obs}.dat')
        for file in sorted(files):
            Path_file = Path(file)
            star_name, analysis = mypy.split_line(Path_file.stem, '_')
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

Path(f'{config.n_sigma_spectrobox}sigmaSpectro_output_tables/').mkdir(parents=True, exist_ok=True)
df_AICc_MD.to_csv(f'{config.n_sigma_spectrobox}sigmaSpectro_output_tables/AICc_values_MD.tsv', sep='\t',index=False)
df_AICc_Chi2.to_csv(f'{config.n_sigma_spectrobox}sigmaSpectro_output_tables/AICc_values_Chi2.tsv', sep='\t',index=False)

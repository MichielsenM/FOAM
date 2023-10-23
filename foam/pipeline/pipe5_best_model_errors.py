"""Calculate the 2 sigma uncertainty region of the maximum likelihood solution using Bayes' theorem."""

import glob, sys
import pandas as pd
import numpy as np
from pathlib import Path
from foam import support_functions as sf
from foam.pipeline.pipeline_config import config
################################################################################
N_dict = config.N_dict  # number of observables
sigma = 2
percentile = {1 : 0.68, 2:0.95, 3:0.997}

if config.observable_additional is not None:
    extra_obs = '+extra'
else:
    extra_obs = ''
################################################################################
def likelihood_chi2(chi2):
    """ Likelihood function of reduced chi-squared """
    return np.exp(-0.5*chi2/(N_dict[obs]-config.k))

def likelihood_MD(MD):
    """ Likelihood function of the mahalanobis distance """
    df_AICc_MD = pd.read_table(f'V_matrix/{config.star}_determinant_conditionNr.tsv', delim_whitespace=True, header=0)
    lndetV = float( (df_AICc_MD.loc[df_AICc_MD['method'] == f'{config.star}_{analysis}', 'ln(det(V))']).iloc[0] )
    return np.exp(-0.5*( MD + config.k*np.log(2*np.pi) + lndetV ))

################################################################################
if config.n_sigma_box != None:
    directory_prefix = f'{config.n_sigma_box}sigmaBox_'
else:
    directory_prefix = f''

for merit in config.merit_functions:
    for obs in config.observable_seismic:
        obs+=extra_obs
        files = glob.glob(f'{directory_prefix}meritvalues/*{merit}_{obs}.hdf')
        for file in sorted(files):
            Path_file = Path(file)
            output_name = Path_file.with_stem(f'{Path_file.stem}_{sigma}sigma-error-ellipse')
            if output_name.is_file(): # Don't duplicate if file is already present
                config.logger.warning(f'file already existed: {output_name}')
                continue

            star_name, analysis = sf.split_line(Path_file.stem, '_')
            df = pd.read_hdf(file)
            df = df.sort_values('meritValue', ascending=True)

            # Dictionary containing different likelihood functions
            switcher={'CS': likelihood_chi2,
                      'MD': likelihood_MD }
            # get the desired function from the dictionary. Returns the lambda function if option is not in the dictionary.
            likelihood_function = switcher.get(merit, lambda x : sys.exit(config.logger.error(f'invalid type of maximum likelihood estimator:{merit}')))

            probabilities = {}
            for column_name in config.free_parameters:
                probabilities.update({column_name:{}})
                for value in df[column_name].unique():     # construct dictionary
                    probabilities[column_name].update({value:0})
                for value in df[column_name]:              # sum over all occurences of parameter values
                    probabilities[column_name][value]+=1
                for value in df[column_name].unique():     # divide by total number of models to get probabilities
                    probabilities[column_name][value] = probabilities[column_name][value] / len(df)

            total_probability=0
            for i in range(len(df)):    # calculate the denominator
                prob = likelihood_function( df.iloc[i]['meritValue']-df.iloc[0]['meritValue'] )
                # prob = likelihood_function( df.iloc[i]['meritValue'] )

                for column_name in config.free_parameters:
                    value = df.iloc[i][column_name]
                    prob = prob*probabilities[column_name][value]
                total_probability+=prob
            p = 0
            for i in range(len(df)):
                prob = likelihood_function( df.iloc[i]['meritValue']-df.iloc[0]['meritValue'] )
                # prob = likelihood_function( df.iloc[i]['meritValue'] )
                for column_name in config.free_parameters:
                    value = df.iloc[i][column_name]
                    prob = prob*probabilities[column_name][value]

                p +=prob/total_probability
                if p>=percentile[sigma]:
                    # Write all models enclosed within the error ellips to a separate file
                    df.iloc[:i+1].to_hdf( output_name , 'models_in_2sigma_error_ellipse', format='table', mode='w')
                    config.logger.debug(f'---------- {analysis} ---------- {i+1} --- {p}')
                    break

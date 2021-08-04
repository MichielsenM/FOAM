"""Calculate the 2 sigma uncertainty region of the maximum likelihood solution.
The figures that are made are not generalised for new grids,
so labels and axis ranges in the beginning of the script will need to be adjusted."""

import glob, sys
import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pypulse import my_python_functions as mypy
import config # imports the config file relative to the location of the main script

################################################################################
N_dict = config.N_dict  # number of observables
merit_abbrev = {'chi2': 'CS', 'mahalanobis': 'MD'}

xlabels_dict = {'M': r'M$_{\rm ini}$', 'Z': r'Z$_{\rm ini}$', 'logD':r'log(D$_{\rm env}$)', 'aov':r'$\alpha_{\rm CBM}$','fov':r'f$_{\rm CBM}$','Xc':r'$\rm X_c$'}
axis_range_dict =  {'M': [2.78, 3.72], 'Z': [0.0147, 0.0233], 'logD':[-0.08, 2.08], 'aov':[-0.011, 0.311],'fov':[-0.0011, 0.0311],'Xc':[0.295, 0.605]}
grid_stepsize_dict =  {'M': 0.1, 'Z': 0.004, 'logD':0.5, 'aov':0.05, 'fov':0.005, 'Xc':0.02}

color_dict = {0:'blue', -1:'blue', -2:'blue', -3:'red', -4:'red', -5:'red'}
ylabel_dict = {0:f'Longest sequence', -1:'Highest amplitude', -2:'Highest frequency', -3:'Longest sequence', -4:'Highest amplitude', -5:'Highest frequency'}

# add ticks and minor ticks so that each step in the grid has at least a minor tick
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
xticks_dict = {'M':[2.8, 3.0, 3.2, 3.4, 3.6], 'Z':[0.015, 0.019, 0.023], 'aov':[0, 0.1, 0.2, 0.3], 'fov':[0, 0.01, 0.02, 0.03], 'logD':[0, 1, 2], 'Xc':[0.3, 0.4, 0.5, 0.6]}
minorticks_dict = {'M':2, 'Z':1, 'aov':2, 'fov':2, 'logD':2, 'Xc':5}

sigma = 2
percentile = {1 : 0.68, 2:0.95, 3:0.997}
################################################################################
def likelihood_chi2(chi2):
    """ Likelihood function of reduced chi-squared """
    return np.exp(-0.5*chi2/(N_dict[obs]-config.k))

def likelihood_MD(MD):
    """ Likelihood function of the mahalanobis distance """
    df_AICc_MD = pd.read_table('V_matrix/determinant_conditionNr.tsv', delim_whitespace=True, header=0)
    lndetV = float( df_AICc_MD.loc[df_AICc_MD['method'] == f'{config.star}_{analysis}', 'ln(det(V))'] )
    return np.exp(-0.5*( MD + config.k*np.log(2*np.pi) + lndetV ))

################################################################################
with open(f'{config.n_sigma_spectrobox}sigmaSpectro_output_tables/{sigma}sigma_errorMargins.txt', 'w') as outfile:
    outfile.write(f'{config.free_param}'+'\n')

for merit in config.merit_functions:
    for obs in config.observable_aic:
        # Make the figure
        width = 4+(2*config.k)
        fig=plt.figure(figsize=(width,2.5))
        gs=GridSpec(ncols=config.k, nrows=1, figure=fig)
        f_ax = {}
        for i in range(config.k):
            if i !=0:
                f_ax.update({ i : fig.add_subplot(gs[0, i], sharey=f_ax[0])})
                f_ax[i].tick_params(axis='y', labelbottom=False, bottom=False, top=False)
                plt.setp(f_ax[i].get_yticklabels(), visible=False)
            else:
                f_ax.update({ i : fig.add_subplot(gs[0, i])})
                f_ax[i].set_yticks([0, -1, -2, -3, -4, -5])

            f_ax[i].set_xlabel(xlabels_dict[config.free_param[i]], size=17)
            f_ax[i].set_xlim(axis_range_dict[config.free_param[i]])
            f_ax[i].tick_params(labelsize=14, length=6)
            f_ax[i].tick_params(which='minor', length=4)

        files = glob.glob(f'{config.n_sigma_spectrobox}sigmaSpectro_extracted_freqs/*{merit_abbrev[merit]}_{obs}.dat')
        j=0
        for file in sorted(files):
            Path_file = Path(file)
            star_name, analysis = mypy.split_line(Path_file.stem, '_')
            df = pd.read_csv(file, delim_whitespace=True, header=0)
            df = df.sort_values('meritValue', ascending=True)

            # Dictionary containing different likelihood functions
            switcher={ 'chi2': likelihood_chi2,
                       'mahalanobis' : likelihood_MD }
            # get the desired function from the dictionary. Returns the lambda function if option is not in the dictionary.
            likelihood_function = switcher.get(merit, lambda x : sys.exit(config.logger.error(f'invalid type of maximum likelihood estimator:{merit}')))

            probabilities = {}
            error_region = {}

            for column_name in config.free_param:
                probabilities.update({column_name:{}})
                error_region.update({column_name:[]})
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

                for column_name in config.free_param:
                    value = df.iloc[i][column_name]
                    prob = prob*probabilities[column_name][value]
                total_probability+=prob
            p = 0
            for i in range(len(df)):
                prob = likelihood_function( df.iloc[i]['meritValue']-df.iloc[0]['meritValue'] )
                # prob = likelihood_function( df.iloc[i]['meritValue'] )
                for column_name in config.free_param:
                    value = df.iloc[i][column_name]
                    prob = prob*probabilities[column_name][value]

                    if value not in error_region[column_name]:
                        error_region[column_name].append(value)
                        error_region[column_name] = sorted(error_region[column_name])

                p +=prob/total_probability
                if p>=percentile[sigma]:
                    config.logger.info(f'---------- {analysis} ---------- {i+1} --- {p}')
                    break

            config.logger.info(error_region)

            with open(f'{config.n_sigma_spectrobox}sigmaSpectro_output_tables/{sigma}sigma_errorMargins.txt', 'a') as outfile:
                outfile.write(f'{analysis} ')
                i=0
                for column_name in config.free_param:
                    outfile.write(f'{error_region[column_name]} ')
                    f_ax[i].scatter(df.iloc[0][column_name], j, color=color_dict[j])
                    f_ax[i].hlines(y=j, xmin=min(error_region[column_name]), xmax=max(error_region[column_name]), color='black')

                    if len(error_region[column_name]) == 1:
                        f_ax[i].hlines(y=j, xmin=max(error_region[column_name][0]-grid_stepsize_dict[column_name], min(xticks_dict[column_name])), xmax=min(error_region[column_name][0]+grid_stepsize_dict[column_name], max(xticks_dict[column_name])), color='grey', alpha=0.7, ls='--')

                    f_ax[i].set_xticks(xticks_dict[column_name])
                    f_ax[i].xaxis.set_minor_locator(AutoMinorLocator(minorticks_dict[column_name]))
                    i+=1
                outfile.write(f'  \\\\ \n')
            j-=1

        labels = [item.get_text() for item in f_ax[0].get_yticklabels()]
        for j in range(len(labels)):
            labels[j] = ylabel_dict[-j]
        f_ax[0].set_yticklabels(labels)

        output_folder = f'figures_errors_{config.n_sigma_spectrobox}sigmaSpectro'
        Path(output_folder).mkdir(parents=True, exist_ok=True)
        leftspace = 0.22-0.015*config.k
        fig.subplots_adjust(left=leftspace, right=0.98, bottom=0.25, top=0.98, wspace=0.2)
        fig.savefig(f'{output_folder}/errors_{merit}_{obs}.png', dpi=300)

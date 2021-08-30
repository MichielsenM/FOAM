""" Functions to perform different kinds of maximum likelihood estimation for the models in a grid, and make correlation plots.
Note: The file with observations needs to hold temperature as Teff, although the analysis is done using the logTeff values."""
# from foam import maximum_likelihood_estimator as mle
import numpy as np
import pandas as pd
import sys, logging, os
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from pathlib import Path
from foam import support_functions as sf
from foam import functions_for_gyre as ffg

logger = logging.getLogger('logger.mle_estimator')  # Make a child logger of "logger" made in the top level script
################################################################################
def plot_correlations(merit_values_file, observations_file, fig_title=None, label_size=20, fig_outputDir='figures_correlation/',
                      percentile_to_show=0.5, logg_or_logL='logL', mark_best_model= False):
    """
    Make a plot of all variables vs each other variable, showing the MLE values as colorscale.
    A kiel/HR diagram is made, depending on if logg_obs or logL_obs is passed as a parameter.
    The subplots on the diagonal show the distribution of that variable.
    The list of variables is retrieved from columns of the merit_values_file,
    where the first column is 'meritValue', which are the MLE values.
    The resulting figure is saved afterwards in the specified location.
    ------- Parameters -------
    merit_values_file: string
        Path to the tsv files with the merit function values and parameters of the models in the grid.
    observations_file: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    fig_title: string
        Title of the figure and name of the saved png.
    label_size: int
        Size of the axis labels.
    fig_outputDir: string
        Output directory for the figures.
    percentile_to_show: float
        Percentile of models to show in the plots.
    logg_or_logL: string
        String 'logg' or 'logL' indicating wheter log of surface gravity (g) or luminosity (L) is plot.
    mark_best_model: boolean
        Indicate the best model with a marker
    """
    # theoretical models
    df_Theo = pd.read_csv(merit_values_file, delim_whitespace=True, header=0)
    df_Theo = df_Theo.sort_values('meritValue', ascending=False)    # Order from high to low, to plot lowest values last
    df_Theo = df_Theo.iloc[int(df_Theo.shape[0]*(1-percentile_to_show)):] # only plot the given percentage lowest meritValues

    df = df_Theo.drop(columns=['rot', 'logTeff', 'logL', 'logg'], errors='ignore') # make new dataframe without the spectroscopic info

    ax_dict={}  # dictionary of dictionaries, holding the subplots of the figure, keys indicate position (row, column) of the subplot
    nr_params = len(df.columns)-1
    for i in range(nr_params):
        ax_dict.update({i:{}})

    fig=plt.figure(figsize=(10,8))
    gs=GridSpec(nr_params,nr_params) # multiple rows and columns
    # proper label format on figures
    axis_labels_dict = {'M': r'M$_{\rm ini}$', 'Z': r'Z$_{\rm ini}$', 'logD':r'log(D$_{\rm env}$)', 'aov':r'$\alpha_{\rm CBM}$','fov':r'f$_{\rm CBM}$','Xc':r'$\rm X_c$'}

    if mark_best_model: min_index = df['meritValue'].idxmin(axis='index', skipna=True)    # get the best model according to the point estimator

    for ix in range(0, nr_params):
        for iy in range(0, nr_params-ix):
            if iy==0:
                shX = None
            else:
                shX = ax_dict[0][ix]
            if (ix==0) or (iy+ix == nr_params-1):
                shY = None
            else:
                shY = ax_dict[iy][0]

            # create subplots and add them to the dictionary
            ax = fig.add_subplot(gs[nr_params-iy-1:nr_params-iy,ix:ix+1], sharex=shX, sharey=shY)
            ax_dict[iy].update({ix:ax})

            # manage visibility and size of the labels and ticks
            ax.tick_params(labelsize=label_size-4)
            if ix == 0:
                ax.set_ylabel(axis_labels_dict[df.columns[iy+1]], size=label_size)
                if (iy == nr_params-1):
                    ax.set_yticklabels(['']*10) # remove ticklabels for just the histogram
            else:
                plt.setp(ax.get_yticklabels(), visible=False)
            if iy == 0:
                ax.set_xlabel(axis_labels_dict[df.columns[nr_params-ix]], size=label_size)
                ax.tick_params(axis='x', rotation=45)
            else:
                plt.setp(ax.get_xticklabels(), visible=False)

            if (iy+ix == nr_params-1):  # make distribution plots on the diagonal subplots
                values = sorted(np.unique(df.iloc[:,nr_params-ix]))
                # determine edges of the bins for the histogram distribution plots
                if len(values) > 1:
                    bin_half_width = (values[0]+values[1])/2-values[0]
                else:
                    bin_half_width = 1E-3
                bin_edges = [values[0]-bin_half_width]
                for i in range(len(values)-1):
                    bin_edges.extend([(values[i]+values[i+1])/2])
                bin_edges.extend([values[-1]+bin_half_width])

                ax.hist( df.iloc[:,nr_params-ix], bins=bin_edges, density=True, cumulative=False, histtype='step' )
                ax.tick_params(axis='y',left=False)
                continue

            im = ax.scatter(df.iloc[:,nr_params-ix], df.iloc[:,iy+1], c=np.log10(df.iloc[:,0]), cmap='hot')
            if mark_best_model: ax.scatter(df.loc[min_index][nr_params-ix], df.loc[min_index][iy+1], color='white', marker = 'x')
            # Adjust x an y limits of subplots
            limit_adjust = (max(df.iloc[:,iy+1]) - min(df.iloc[:,iy+1]))*0.08
            ax.set_ylim( min(df.iloc[:,iy+1])-limit_adjust,  max(df.iloc[:,iy+1])+limit_adjust  )
            limit_adjust = (max(df.iloc[:,nr_params-ix]) - min(df.iloc[:,nr_params-ix])) *0.08
            ax.set_xlim( min(df.iloc[:,nr_params-ix])-limit_adjust, max(df.iloc[:,nr_params-ix])+limit_adjust )

    fig.align_labels()

    # add subplot in top right for Kiel or HRD
    ax_hrd = fig.add_axes([0.508, 0.65, 0.33, 0.33]) # X, Y, widht, height

    ax_hrd.set_xlabel(r'log(T$_{\mathrm{eff}}$ [K])', size=label_size)
    ax_hrd.tick_params(labelsize=label_size-4)
    ax_hrd.invert_xaxis()

    im = ax_hrd.scatter(df_Theo['logTeff'], df_Theo[logg_or_logL], c=np.log10(df_Theo['meritValue']), cmap='hot')
    ax_hrd.set_ylabel(f'{logg_or_logL[:-1]} {logg_or_logL[-1]}')
    if logg_or_logL == 'logL':
        ax_hrd.set_ylabel(r'log(L [L$_{\odot}$])', size=label_size)
    elif logg_or_logL == 'logg':
        ax_hrd.set_ylabel(r'log$g$ [dex]', size=label_size)

    # observations
    Obs_dFrame  = pd.read_table(observations_file, delim_whitespace=True, header=0)
    # Observed spectroscopic error bar
    # To add the 1 and 3 sigma spectro error boxes, calculate their width (so 2 and 6 sigmas wide)
    width_logTeff_2sigma= np.log10(Obs_dFrame['Teff'][0]+Obs_dFrame['Teff_err'][0]) - np.log10(Obs_dFrame['Teff'][0]-Obs_dFrame['Teff_err'][0])
    width_logTeff_6sigma= np.log10(Obs_dFrame['Teff'][0]+3*Obs_dFrame['Teff_err'][0]) - np.log10(Obs_dFrame['Teff'][0]-3*Obs_dFrame['Teff_err'][0])
    errorbox_1s = patches.Rectangle((np.log10(Obs_dFrame['Teff'][0]-Obs_dFrame['Teff_err'][0]),Obs_dFrame[logg_or_logL][0]-Obs_dFrame[f'{logg_or_logL}_err'][0]), width_logTeff_2sigma, 2*Obs_dFrame[f'{logg_or_logL}_err'][0],linewidth=1.7,edgecolor='cyan',facecolor='none')
    errorbox_3s = patches.Rectangle((np.log10(Obs_dFrame['Teff'][0]-3*Obs_dFrame['Teff_err'][0]),Obs_dFrame[logg_or_logL][0]-3*Obs_dFrame[f'{logg_or_logL}_err'][0]), width_logTeff_6sigma, 6*Obs_dFrame[f'{logg_or_logL}_err'][0],linewidth=1.7,edgecolor='cyan',facecolor='none')
    ax_hrd.add_patch(errorbox_1s)
    ax_hrd.add_patch(errorbox_3s)
    if mark_best_model: ax_hrd.scatter(df_Theo['logTeff'][min_index], df_Theo[logg_or_logL][min_index], marker='x', color='white')

    # Add color bar
    cax = fig.add_axes([0.841, 0.266, 0.05, 0.6]) # X, Y, widht, height
    cbar= fig.colorbar(im, cax=cax, orientation='vertical')
    if '_MD_' in fig_title:
        cbar.set_label('log(MD)', rotation=90, size=label_size)
    elif '_CS_' in fig_title:
        cbar.set_label(r'log($\chi^2$)', rotation=90, size=label_size)
    else:
        cbar.set_label('log(merit function value)', rotation=90)
    cbar.ax.tick_params(labelsize=label_size-4)
    fig.subplots_adjust(left=0.109, right=0.835, bottom=0.13, top=0.99)

    # fig.suptitle(fig_title, horizontalalignment='left', size=20, x=0.28)
    Path(fig_outputDir).mkdir(parents=True, exist_ok=True)
    fig.savefig(f'{fig_outputDir}{fig_title}.png', dpi=400)
    plt.close('all')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def calculate_likelihood(Obs_path, Theo_file, observables=None, merit_function =None, star_name=None):
    """
    Perform a maximum likelihood estimation using the provided type of merit function on the list of  observables.
    Writes a data file with the values of the merit funtion and input parameters of each model.
    ------- Parameters -------
    Obs_path: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    Theo_file: string
        Path to the tsv file with the theoretical model input parameters (first set of columns), frequency or period values (last set of columns),
        and possibly extra columns with additional observables (these columns should be in between the input parameters and frequency columns).
    observables: list of strings
        Can contain 'frequencies' or 'periods', 'period_spacing', and 'rope_length', which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.
    merit_function: string
        The type of merit function to use. Currently supports "chi2" and "mahalanobis".
    star_name: string
        Name of the star, used in file naming.
    """
    # Read in the observed data and make an array of the observed obervables
    Obs_dFrame = pd.read_table(Obs_path, delim_whitespace=True, header=0)
    Obs, ObsErr, file_suffix_observables = create_obs_observables_array(Obs_dFrame, observables)

    Path_theo   = Path(Theo_file)
    #suffix for filename to indicate the merit funtion used
    suffix = {'chi2'       : 'CS',
              'mahalanobis': 'MD'}

    # set the name of the output file
    head, tail = sf.split_line(Path_theo.stem, star_name)
    DataOut = f'{Path_theo.parent}/{star_name}{tail}_{suffix[merit_function]}_{file_suffix_observables}.dat'

    # Theoretical grid data
    Theo_dFrame = pd.read_table(Theo_file, delim_whitespace=True, header=0)
    Thetas      = np.asarray(Theo_dFrame.loc[:,:'Xc']) # varied parameters in the grid (e.g. Mini, Xini, Xc etc.)
    Theo_puls   = np.asarray(Theo_dFrame.loc[:,'f1':]) # theoretical pulsations corresponding to the observed ones

    missing_absolute = np.where(Theo_dFrame.columns.to_series().str.contains('f_missing'))[0]               # get the interruptions in the pattern, absolute index in dataframe
    missing_relative = np.where(Theo_dFrame.loc[:,'f1':].columns.to_series().str.contains('f_missing'))[0]  # get the interruptions in the pattern, index relative within pulsations
    Theo_dFrame = Theo_dFrame.drop(columns=Theo_dFrame.columns[missing_absolute])   # Remove columns of missing frequencies
    missing_indices=[ missing_relative[i]-i for i in range(len(missing_relative)) ]         # Adjust indices for removed lines of missing frequencies

    # Make new list of theoretical models without entries with value -1
    newTheo = []
    newThetas = []
    for i in range(len(Theo_puls)):
        if (Theo_puls[i][0] !=-1) and (Theo_puls[i][-1] !=-1):     # ignore models where one of the freqs is -1
            theo_observables = create_theo_observables_array(Theo_dFrame, i, observables, missing_indices)   # make an array of the theoretical observables for each model

            newTheo.append(theo_observables)
            newThetas.append(Thetas[i])
    neg_value = np.unique(np.where(Theo_puls==-1)[0])
    logger.debug(f'File: {Theo_file}')
    logger.debug(f'Observables      : {observables}')
    logger.debug(f'ignored -1 freqs : {len(neg_value)}')
    logger.debug(f'original #models : {len(Theo_puls)}')
    logger.debug(f'remaining #models: {len(newTheo)}')

    Theo_observables = np.asarray(newTheo)
    Thetas           = np.asarray(newThetas)

    # Dictionary containing different merit functions
    switcher={ 'chi2': merit_chi2,
                'mahalanobis' : merit_mahalanobis}

    # get the desired function from the dictionary. Returns the lambda function if option is not in the dictionary.
    selected_merit_function = switcher.get(merit_function, lambda x, y, z: sys.exit(logger.error('invalid type of maximum likelihood estimator')))
    merit_values = selected_merit_function(Obs, ObsErr, Theo_observables, fig_title=f'{star_name}{tail}_{suffix[merit_function]}_{file_suffix_observables}')

    # Print smallest and highest values
    idx2 = np.argsort(merit_values)
    Parameters = Theo_dFrame.loc[:,:'Xc'].columns  # Parameter names
    logger.info(f'Smallest {merit_function} : {merit_values[np.argsort(merit_values)][0]}')
    logger.info(f'Highest {merit_function}  : {merit_values[np.argsort(merit_values)][-1]}')
    logger.info(f'meritValue & {Parameters[0]}   & {Parameters[1]}  & {Parameters[2]}  &{Parameters[3]}& {Parameters[4]} &  {Parameters[5]} & {Parameters[6]}')
    logger.info('-------------------------------------------------------')
    # Print the ten models with the smallest values
    for i in range(10):
        row = f'{merit_values[np.argsort(merit_values)][i]:.4f}'
        for k in range(np.shape(Thetas[idx2,:])[1]):
            row += f' & {Thetas[idx2,:][i,k]:.3f}'
        logger.info(row)

    # Save the results
    CombData = np.concatenate((np.matrix(merit_values).T,Thetas),axis=1)  # add an additional column for MLE 'meritValues'
    Parameters= Parameters.insert(0, 'meritValue') # add an additional parameter name

    df = pd.DataFrame(data=CombData, columns=Parameters) # put the data in a pandas DataFrame

    for spectro in ['logTeff', 'logL', 'logg']:
        if spectro in Theo_dFrame.columns:
            df = pd.merge(df, Theo_dFrame[['Z', 'M', 'logD', 'aov', 'fov', 'Xc', spectro]], how='inner', on=['Z', 'M', 'logD', 'aov', 'fov', 'Xc'])
    df.to_csv(DataOut, sep='\t',index=False)             # write the dataframe to a tsv file
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def create_theo_observables_array(Theo_dFrame, index, observables_in, missing_indices):
    """
    Create an array of theoretical observables.
    ------- Parameters -------
    Theo_dFrame: pandas dataFrame
        DataFrame containing the theoretical periods or frequencies (as the last columns), along with any additional
        columns containing extra observables.
    index: int
        Row index in the dataFrame of the theoretical model to make the array for.
    observables_in: list of strings
        Which observables are included in the returned array.
        Can contain 'frequency', 'period', 'period_spacing', and/or 'rope_length', which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.
    missing_indices: list of int
        Contains the indices of the missing pulsations so that the period sapcing pattern can be split around them.
    ------- Returns -------
    observables_out: numpy array of floats
        The values of the specified observables for the model.
    """
    observables = list(observables_in)  # Make a copy to leave the array handed to this function unaltered.
    observables_out = np.asarray(Theo_dFrame.loc[index,'f1':])  # add the periods or frequencies to the output list

    if 'period' in observables:
        periods = np.asarray(Theo_dFrame.loc[index,'f1':])      # a separate list of periods that is preserved after adding other observables
        observables.remove('period')

    elif 'frequency' in observables:
        periods = 1/np.asarray(Theo_dFrame.loc[index,'f1':])      # a separate list of periods that is preserved after adding other observables
        observables.remove('frequency')
    else:
        periods = np.asarray(Theo_dFrame.loc[index,'f1':])  # Assume the tsv file was in periods if nothing was specified
        observables_out = np.asarray([])                    # Don't use period or freq as observables, so overwrite previous list to be empty

    if 'period_spacing' in observables:
        for periods_part in np.split(periods,missing_indices):
            spacing = ffg.generate_thry_series(periods_part)
            spacing = np.asarray(spacing)/86400 # switch back from seconds to days (so both P and dP are in days)
            observables_out = np.append(observables_out, spacing)   # Include dP as observables
        observables.remove('period_spacing')

    if 'rope_length' in observables:
        for periods_part in np.split(periods,missing_indices):
            RL, error = PdP_pattern_rope_length(periods_part)
            observables_out = np.append(observables_out, RL)        # Include pattern rope length as observable
        observables.remove('rope_length')

    # Add all other observables in the list from the dataFrame
    for observable in observables:
        observables_out = np.append(observables_out, np.asarray(Theo_dFrame.loc[index,observable]))

    return observables_out

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def create_obs_observables_array(Obs_dFrame, observables):
    """
    Create an array of the observed observables.
    ------- Parameters -------
    Obs_dFrame: pandas dataFrame
        DataFrame containing the theoretical frequencies, periods, and any additional observables as columns, as well as columns with their errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    observables: list of strings
        Which observables are included in the returned array, must contain 'frequency' or 'period'.
        Can contain 'period_spacing' and 'rope_length', which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.

    ------- Returns -------
    observables_out, observablesErr_out: numpy array of floats
        The values of the specified observables (or their errors).
    filename_suffix: string
        suffix for the filename, containing all the included observables, separated by '-'
    """
    missing_indices = np.where(Obs_dFrame.index.isin(['f_missing']))[0] # get the interruptions in the pattern
    missing_indices=[ missing_indices[i]-i for i in range(len(missing_indices)) ]       # Adjust indices for removed lines of missing frequencies
    if len(missing_indices)!=0: Obs_dFrame = Obs_dFrame.drop(index='f_missing')  # remove lines indicating missing frequencies (if they are present)

    observables=list(observables)  #make a copy of the list, to not alter the one that was given to the function
    period = np.asarray(Obs_dFrame['period'])
    periodErr = np.asarray(Obs_dFrame['period_err'])
    filename_suffix = ''

    periods_parts = np.split(period,missing_indices)
    periodsErr_parts = np.split(periodErr,missing_indices)

    if 'period' in observables:
        observables_out = np.asarray(Obs_dFrame['period'])
        observablesErr_out = np.asarray(Obs_dFrame['period_err'])
        filename_suffix = 'P'
        observables.remove('period')

    elif 'frequency' in observables:
        observables_out = np.asarray(Obs_dFrame['frequency'])
        observablesErr_out = np.asarray(Obs_dFrame['frequency_err'])
        filename_suffix = 'f'
        observables.remove('frequency')
    else:
        observables_out = np.asarray([])
        observablesErr_out = np.asarray([])

    if 'period_spacing' in observables:
        for periods, periodsErr in zip(periods_parts, periodsErr_parts):
            spacing, spacing_errs = ffg.generate_obs_series(periods, periodsErr)
            spacing = np.asarray(spacing)/86400 # switch back from seconds to days (so both P and dP are in days)
            spacing_errs = np.asarray(spacing_errs)/86400

            observables_out = np.append(observables_out, spacing)   # Include dP as observables
            observablesErr_out = np.append(observablesErr_out, spacing_errs)

        if filename_suffix != '': filename_suffix+='-'     # only add - if dP is not first observable
        filename_suffix+='dP'
        observables.remove('period_spacing')

    if 'rope_length' in observables:
        for periods, periodsErr in zip(periods_parts, periodsErr_parts):
            RL, error = PdP_pattern_rope_length(periods, P_error=periodsErr)
            observables_out = np.append(observables_out, RL)      # Include pattern rope length as observable
            observablesErr_out = np.append(observablesErr_out, error)

        filename_suffix+='-rl'
        observables.remove('rope_length')

    # Add all other observables in the list from the dataFrame
    for observable in observables:
        if observable == 'logTeff': observable = 'Teff' # To read it as Teff from the observations datafile
        Obs = np.asarray(Obs_dFrame[observable]) # since these columns have less entries, the dataFrame has NaNs in the empty rows
        Obs = Obs[~np.isnan(Obs)]                # remove all entries that are NaN in the numpy array
        ObsErr = np.asarray(Obs_dFrame[f'{observable}_err'])
        ObsErr = ObsErr[~np.isnan(ObsErr)]

        if observable == 'Teff': observable = 'logTeff' # To write it as logTeff in the filename
        observables_out = np.append(observables_out, Obs)
        observablesErr_out = np.append(observablesErr_out, ObsErr)
        filename_suffix+=f'-{observable}'

    return observables_out, observablesErr_out, filename_suffix

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def merit_chi2(YObs, ObsErr, YTheo, fig_title=None):
    """
    Calculate chi squared values for the given theoretial patterns
    ------- Parameters -------
    YObs, ObsErr: numpy array of floats
        Observed values and their errors (period or frequency)
    YTheo: numpy array of arrays of floats
        Array of all theoretical patterns to calculate the chi squared value for.

    ------- Returns -------
    chi2: numpy array of floats
        chi squared values for the given theoretical values
    """
    chi2 = np.array([  np.sum(((one_YTheo-YObs)/ObsErr)**2) for one_YTheo in YTheo ])
    return chi2

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def merit_mahalanobis(YObs, ObsErr, YTheo, fig_title=None):
    """
    Calculate mahalanobis distance values for the given theoretial patterns
    ------- Parameters -------
    YObs, ObsErr: numpy array of floats
        Observed values and their errors (period or frequency)
    YTheo: numpy array of arrays of floats
        Array of all theoretical patterns to calculate the chi squared value for.

    ------- Returns -------
    MD: numpy array of floats
        Mahalanobis distances for the given theoretical values
    """
    # Convert to matrix format
    YObsMat = np.matrix(YObs).T
    YTheoMat = np.matrix(YTheo).T

    # Calculate the average on the theoretical values (e.g. frequencies)
    # over the entire grid. Returns a vector of dimension N x 1
    Yav = YTheoMat.mean(1)

    # Calculate the variance-covriance matrix
    q = np.shape(YTheo)[0]  # number of grid points
    N = len(YObs)           # number of observed values
    V = np.zeros((N,N))
    for i in range(q):
        difference = np.subtract(YTheoMat[:,i],Yav)
        V += np.matmul(difference,difference.T)
    V = V/float(q-1)

    # Include observational errors in the variance-covariance matrix
    V = V + np.diag(ObsErr**2.)
    check_matrix(V, fig_title=fig_title)     # check if positive definite and make figure
    # Calculate Mahalanobis distances
    MD = np.zeros(q)
    Vinv = np.linalg.inv(V)
    for i in range(q):
        diff = (YTheoMat[:,i]-YObsMat)
        MD[i] = np.matmul(np.matmul(diff.T,Vinv),diff)[0][0]

    return MD

################################################################################
def check_matrix(V, plot=True, fig_title='Vmatrix'):
    """
    Check the if the the eigenvalues of the Variance-covariance matrix are all positive,
    since this means the matrix is positive definite. Compute its determinant and condition number,
    and write them to a tsv file. Create and save a figure of the variance-covariance matrix.
    ------- Parameters -------
    V: 2D np array
        Variance-covariance matrix
    plot: boolean
        Flag to make a plot of the variance-covariance matrix
    fig_title: string
        The name of the figure to be created.
    """
    if np.all(np.linalg.eigvals(V) > 0)==False: # If all eigencalues are >0, it is positive definite
        sys.exit(logger.error('V matrix is possibly not positive definite (since eigenvalues are not all > 0)'))

    logger.info(f'max(V) = {np.max(V)}')
    kk=10 # multiply the matrix by the exponent of this, otherwise the determinant can be too small for the numerics
    file_Path = Path(f'{os.getcwd()}/V_matrix/determinant_conditionNr.tsv')
    file_Path.parent.mkdir(parents=True, exist_ok=True)

    if not file_Path.is_file():
        with file_Path.open("w") as file:
            file.write(f'method \t ln(det(V)) \t condition_number \n')
    with file_Path.open("a") as file:
        file.write(f'{fig_title} \t {np.log(np.linalg.det(np.exp(kk)*V))-kk*V.shape[0]:.2f} \t {np.linalg.cond(V):.2f} \n ')

    if plot is True:
        im = plt.imshow(V*10**4, aspect='auto', cmap='Reds') # Do *10^4 to get rid of small values, and put this in the colorbar label
        plt.ylabel(rf'Obs {V.shape[0]}  $\leftarrow$ Obs 1', size=14)
        plt.xlabel(rf'Obs 1 $\rightarrow$ Obs {V.shape[0]} ', size=14)
        # if (V.shape[0]==36):
        #     plt.ylabel(rf'Mode Period {V.shape[0]}  $\leftarrow$ Mode Period 1', size=14)
        #     plt.xlabel(rf'Mode Period 1 $\rightarrow$ Mode Period {V.shape[0]} ', size=14)
        # else:
        #     plt.ylabel(rf'$\Delta P_{ {V.shape[0]} }  \leftarrow \Delta P_{1}$', size=14)
        #     plt.xlabel(rf'$\Delta P_{1} \rightarrow \Delta P_{ {V.shape[0]} }$', size=14)

        cbar = plt.colorbar(im)
        cbar.ax.set_ylabel(r'[d$^{2} 10^{-4}$]', rotation=90, labelpad=15, size=14)
        cbar.ax.tick_params(labelsize=14)
        plt.title('Variance covariance matrix')
        plt.tick_params(labelsize=14)
        plt.tight_layout()
        plt.savefig(f'{os.getcwd()}/V_matrix/{fig_title}.png')
        # plt.savefig(f'{os.getcwd()}/V_matrix/{fig_title}.pdf')
        plt.close('all')

################################################################################
def PdP_pattern_rope_length(P, P_error=[-1]):
    """
    Calculate the rope length of the period spacing pattern (total eucledian distance between each consecutive point in P-dP space).
    ------- Parameters -------
    P, P_error: numpy array of float
        Periods (P) their errors.

    ------- Returns -------
    total_length, total_error: float
        The sum of the eucledian distance between each consecutive point in P-dP space (P in days, dP in seconds)
        and the error on this total length.
    """
    if P_error[0] != -1:
        dP, dP_error = ffg.generate_obs_series(P, P_error)
    else:
        dP = ffg.generate_thry_series(P)

    total_length=0
    deltas_lengths = []
    for i in range(len(dP)-1):
        dx = P[i]-P[i+1]
        dy = dP[i]-dP[i+1]
        total_length += np.sqrt(dx**2+dy**2)

        # calculate the error on the length between each of the P-dP points
        if P_error[0] != -1:
            delta_dx = np.sqrt( P_error[i]**2 + P_error[i+1]**2  )
            delta_dy = np.sqrt( dP_error[i]**2 + dP_error[i+1]**2  )

            delta_length = np.sqrt( ((dx*delta_dx)**2 + (dy*delta_dy)**2) /(dx**2 + dy**2) )
            deltas_lengths.append(delta_length)

    # calculate the total resulting error on the length of the whole pattern
    total_error = 0
    if P_error[0] != -1:
        for delta in deltas_lengths:
            total_error+= delta**2
        total_error = np.sqrt(total_error)

    return total_length, total_error

################################################################################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
def spectro_cutoff(merit_values_file, observations_file, nsigma=2):
    """
    Make an n-sigma cutoff of the models based on the spectoscopic observations.
    Save this as a file with prefix "clipped".
    ------- Parameters -------
    merit_values_file: string
        Path to the tsv files with the merit funtion values and the spectroscopic info of the models in the grid.
    observations_file: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    nsigma: int
        How many sigmas you want to make the interval to accept models.
    """
    Obs_dFrame = pd.read_table(observations_file, delim_whitespace=True, header=0)
    df_Theo = pd.read_csv(merit_values_file, delim_whitespace=True, header=0)

    df_Theo = df_Theo[df_Theo.logTeff < np.log10(Obs_dFrame['Teff'][0]+nsigma*Obs_dFrame['Teff_err'][0])]
    df_Theo = df_Theo[df_Theo.logTeff > np.log10(Obs_dFrame['Teff'][0]-nsigma*Obs_dFrame['Teff_err'][0])]
    df_Theo = df_Theo[df_Theo.logg < Obs_dFrame['logg'][0]+nsigma*Obs_dFrame['logg_err'][0]]
    df_Theo = df_Theo[df_Theo.logg > Obs_dFrame['logg'][0]-nsigma*Obs_dFrame['logg_err'][0]]
    df_Theo = df_Theo[df_Theo.logL < Obs_dFrame['logL'][0]+nsigma*Obs_dFrame['logL_err'][0]]
    df_Theo = df_Theo[df_Theo.logL > Obs_dFrame['logL'][0]-nsigma*Obs_dFrame['logL_err'][0]]

    outputFile = f'{nsigma}sigmaSpectro_{merit_values_file}'
    Path(outputFile).parent.mkdir(parents=True, exist_ok=True)
    df_Theo.to_csv(outputFile, sep='\t',index=False)

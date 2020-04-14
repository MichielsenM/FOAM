""" Functions to perform different kinds of maximum likelihood estimation for the models in a grid, and make correlation plots."""
import numpy as np
import pandas as pd
import sys, logging
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
from . import my_python_functions as mypy
from . import functions_for_gyre as ffg

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('logger')
logger.setLevel(logging.DEBUG)
################################################################################
def plot_correlations(MLE_values_file, observations_file, grid_spectroscopy, fig_title=None, label_size=18, fig_outputDir='figures/',
                      percentile_to_show=0.1, logg_or_logL='logL'):
    """
    Make a plot of all variables vs each other variable, showing the MLE values as colorscale.
    A kiel/HR diagram is made, depending on if logg_obs or logL_obs is passed as a parameter.
    The subplots on the diagonal show the distribution of that variable.
    The list of variables is retrieved from columns of the MLE_values_file,
    where the first column is 'distance', which are the MLE values.
    The resulting figure is saved afterwards in the specified location.
    ------- Parameters -------
    MLE_values_file, grid_spectroscopy: string
        Path to the tsv files with the MLE values / the spectroscopic info of the models in the grid.
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
    """
    # theoretical models
    df = pd.read_csv(MLE_values_file, delim_whitespace=True, header=0)
    df = df.sort_values('distance', ascending=False)    # Order from high to low, to plot lowest values last
    df = df.iloc[int(df.shape[0]*(1-percentile_to_show)):] # only plot the given percentage lowest distances

    ax_dict={}
    nr_params = len(df.columns)-1
    for i in range(nr_params):
        ax_dict.update({i:{}})

    fig=plt.figure(figsize=(10,8))
    gs=GridSpec(nr_params,nr_params) # multiple rows and columns

    ax_hrd = fig.add_subplot(gs[0:3,nr_params-3:nr_params], sharex=None, sharey=None)   # add subplot in top right for Kiel or HRD
    ax_hrd.set_xlabel('Teff')
    ax_hrd.invert_xaxis()
    # ax_hrd.set_xscale("log")

    # get the spectroscopic info and merge the dataFrames into a new one containing all the info of the models they have in common
    spectro_df = pd.read_csv(grid_spectroscopy, delim_whitespace=True, header=0)
    df_merged = pd.merge(df, spectro_df, how='inner', on=['Z', 'M', 'logD', 'aov', 'fov', 'Xc'])

    # observations
    Obs_dFrame  = pd.read_table(observations_file, delim_whitespace=True, header=0)

     # Observed pectroscopic error bar
    ax_hrd.errorbar(Obs_dFrame['Teff'][0], Obs_dFrame[logg_or_logL][0], xerr=Obs_dFrame['Teff_err'][0], yerr=Obs_dFrame[f'{logg_or_logL}_err'][0])
    im = ax_hrd.scatter(df_merged['Teff'], df_merged[logg_or_logL], c=np.log10(df_merged['distance']), cmap='hot')
    ax_hrd.set_ylabel(f'{logg_or_logL[:-1]} {logg_or_logL[-1]}')

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
                ax.set_ylabel(df.columns[iy+1], size=label_size)
                if (iy+ix == nr_params-1):
                    ax.set_yticklabels(['']*10)
            else:
                plt.setp(ax.get_yticklabels(), visible=False)
            if iy == 0:
                ax.set_xlabel(df.columns[nr_params-ix], size=label_size)
            else:
                plt.setp(ax.get_xticklabels(), visible=False)

            if (iy+ix == nr_params-1):  # make distribution plots on the diagonal subplots
                values = sorted(np.unique(df.iloc[:,nr_params-ix]))
                # determine edges of the bins for the histogram distribution plots
                bin_edges = [values[0]-1E-3]
                for i in range(len(values)-1):
                    bin_edges.extend([(values[i]+values[i+1])/2])
                bin_edges.extend([values[-1]+1E-3])

                ax.hist( df.iloc[:,nr_params-ix], bins=bin_edges, density=True, cumulative=False, histtype='step' )
                continue

            im = ax.scatter(df.iloc[:,nr_params-ix], df.iloc[:,iy+1], c=np.log10(df.iloc[:,0]), cmap='hot') # gist_rainbow as alternative cmap
            # Adjust x an y limits of subplots
            ax.set_ylim( max(ax.get_ylim()[0], min(df.iloc[:,iy+1])-0.01 ) ,  min(ax.get_ylim()[1], max(df.iloc[:,iy+1])+0.01 ) )
            ax.set_xlim( max(ax.get_xlim()[0], min(df.iloc[:,nr_params-ix])-0.01 ) ,  min(ax.get_xlim()[1], max(df.iloc[:,nr_params-ix])+0.01 ) )

    fig.align_labels()

    cax = fig.add_axes([0.88, 0.3, 0.05, 0.6]) # X, Y, widht, height
    cbar= fig.colorbar(im, cax=cax, orientation='vertical')
    cbar.set_label('log(MLE value)', rotation=90)
    plt.subplots_adjust(left=0.11, right=0.87, bottom=0.1, top=0.95)

    fig.suptitle(fig_title)
    Path(fig_outputDir).mkdir(parents=True, exist_ok=True)
    fig.savefig(f'{fig_outputDir}{fig_title}.png', dpi=300)
    plt.close('all')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def estimate_max_likelihood(Obs_path, Theo_file, observables=[], estimator_type =''):
    """
    Perform a maximum likelihood estimation using the provided type of estimator on the list of  observables.
    Writes a data file with the MLE values and input parameters of each model.
    ------- Parameters -------
    Obs_path: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    Theo_file: string
        Path to the tsv file with the theoretical model input parameters (first set of columns), frequency or period values (last set of columns),
        and possibly extra columns with additional observables (these columns should be in between the input parameters and frequency columns).
    observables: list of strings
        Which observables are used in the analysis, must contain 'frequencies' or 'periods'.
        Can contain 'period_spacing' and 'rope_length', which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.
    estimator_type: string
        The type of maximum likelihood estimation to use. Currently supports "chi2", "mahalanobis" and "rope_length".
    """
    if estimator_type == 'rope_length' and observables!=['period']:
        sys.exit(logger.error('Can currently only use period as observables when using the rope length as estimator.'))

    # Read in the observed data and make an array of the observed obervables
    Obs_dFrame = pd.read_table(Obs_path, delim_whitespace=True, header=0)
    Obs, ObsErr, file_suffix_observables = create_obs_observables_array(Obs_dFrame, observables)

    Path_theo   = Path(Theo_file)
    #suffix for filename to indicate the MLE method used
    suffix = {'chi2'       : 'CS',
              'mahalanobis': 'MD',
              'rope_length': 'RL'}

    # set the name of the output file
    head, tail = mypy.split_line(Path_theo.stem, 'KIC')
    DataOut = f'{Path_theo.parent}/KIC{tail}_{suffix[estimator_type]}_{file_suffix_observables}.dat'

    # Theoretical grid data
    Theo_dFrame = pd.read_table(Theo_file, delim_whitespace=True, header=0)
    Thetas      = np.asarray(Theo_dFrame.loc[:,:'Xc']) # varied parameters in the grid (e.g. Mini, Xini, Xc etc.)
    Theo_puls   = np.asarray(Theo_dFrame.loc[:,'f1':]) # theoretical pulsations corresponding to the observed ones

    # Make new list of theoretical models without entries with value -1
    newTheo = []
    newThetas = []
    for i in range(len(Theo_puls)):
        if (Theo_puls[i][0] !=-1) and (Theo_puls[i][-1] !=-1):     # ignore models where one of the freqs is -1
            theo_observables = create_theo_observables_array(Theo_dFrame, i, observables)   # make an array of the theoretical observables for each model

            newTheo.append(theo_observables)
            newThetas.append(Thetas[i])
    neg_value = np.unique(np.where(Theo_puls==-1)[0])
    logger.debug(f'ignored -1 freqs : {len(neg_value)}')
    logger.debug(f'original #models : {len(Theo_puls)}')
    logger.debug(f'remaining #models: {len(newTheo)}')

    Theo_observables = np.asarray(newTheo)
    Thetas           = np.asarray(newThetas)

    # Dictionary containing different functions for MLE
    switcher={ 'chi2': mle_chi2,
                'mahalanobis' : mle_mahalanobis,
                'rope_length': mle_rope_length }

    # get the desired function from the dictionary. Returns the lambda function if option is not in the dictionary.
    mle_function = switcher.get(estimator_type, lambda x, y, z: sys.exit(logger.error('invalid type of maximum likelihood estimator')))
    mle_values = mle_function(Obs, ObsErr, Theo_observables)

    # Print smallest and highest values
    idx2 = np.argsort(mle_values)
    Parameters = Theo_dFrame.loc[:,:'Xc'].columns  # Parameter names
    logger.info(f'Smallest {estimator_type} : {mle_values[np.argsort(mle_values)][0]}')
    logger.info(f'Highest {estimator_type}  : {mle_values[np.argsort(mle_values)][-1]}')
    logger.info(f'distance & {Parameters[0]}   & {Parameters[1]}  & {Parameters[2]}  &{Parameters[3]}& {Parameters[4]} &  {Parameters[5]} & {Parameters[6]}')
    logger.info('-------------------------------------------------------')
    # Print the ten models with the smallest values
    for i in range(10):
        row = f'{mle_values[np.argsort(mle_values)][i]:.4f}'
        for k in range(np.shape(Thetas[idx2,:])[1]):
            row += f' & {Thetas[idx2,:][i,k]:.3f}'
        logger.info(row)

    # Save the results
    CombData = np.concatenate((np.matrix(mle_values).T,Thetas),axis=1)  # add an additional column for MLE 'distances'
    Parameters= Parameters.insert(0, 'distance') # add an additional parameter name

    df = pd.DataFrame(data=CombData, columns=Parameters) # put the data in a pandas DataFrame
    df.to_csv(DataOut, sep='\t',index=False)             # write the dataframe to a tsv file

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def create_theo_observables_array(Theo_dFrame, index, observables):
    """
    Create an array of theoretical observables.
    ------- Parameters -------
    Theo_dFrame: pandas dataFrame
        DataFrame containing the theoretical periods or frequencies (as the last columns), along with any additional
        columns containing extra observables.
    index: int
        Row index in the dataFrame of the theoretical model to make the array for.
    observables: list of strings
        Which observables are included in the returned array, must contain 'frequencies' or 'periods'.
        Can contain 'period_spacing' and 'rope_length', which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.

    ------- Returns -------
    observables_out: numpy array of floats
        The values of the specified observables for the model.
    """
    observables=list(observables)  #make a copy of the list, to not alter the one that was given to the function
    observables_out  = np.asarray(Theo_dFrame.loc[index,'f1':])

    if 'period' in observables:
        periods = np.asarray(Theo_dFrame.loc[index,'f1':])      # a separate list of periods that is preserved after adding other observables
        observables.remove('period')

    elif 'frequency' in observables:
        periods = 1/np.asarray(Theo_dFrame.loc[index,'f1':])      # a separate list of periods that is preserved after adding other observables
        observables.remove('frequency')
    else:
        sys.exit(logger.error(f'\"period\" or \"frequency\" should be one of the observables' ))

    if 'period_spacing' in observables:
        spacing = ffg.generate_thry_series(periods)
        spacing = np.asarray(spacing)/86400 # switch back from seconds to days (so both P and dP are in days)
        observables_out = np.append(observables_out, spacing)   # Include dP as observables
        observables.remove('period_spacing')

    if 'rope_length' in observables:
        RL, error = PdP_pattern_rope_length(periods)
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
    observables=list(observables)  #make a copy of the list, to not alter the one that was given to the function
    periods = np.asarray(Obs_dFrame['period'])
    periodsErr = np.asarray(Obs_dFrame['period_err']*4)
    filename_suffix = ''

    if 'period' in observables:
        observables_out = np.asarray(Obs_dFrame['period'])
        observablesErr_out = np.asarray(Obs_dFrame['period_err']*4)   # times 4 to take correlation structure of the data into account (see Moravveji et al. 2016)
        filename_suffix = 'P'
        observables.remove('period')

    elif 'frequency' in observables:
        observables_out = np.asarray(Obs_dFrame['frequency'])
        observablesErr_out = np.asarray(Obs_dFrame['frequency_err']*4)   # times 4 to take correlation structure of the data into account (see Moravveji et al. 2016)
        filename_suffix = 'f'
        observables.remove('frequency')
    else:
        sys.exit(logger.error(f'\"period\" or \"frequency\" should be one of the observables' ))

    if 'period_spacing' in observables:
        spacing, spacing_errs = ffg.generate_obs_series(periods, periodsErr)
        spacing = np.asarray(spacing)/86400 # switch back from seconds to days (so both P and dP are in days)
        spacing_errs = np.asarray(spacing_errs)/86400

        observables_out = np.append(observables_out, spacing)   # Include dP as observables
        observablesErr_out = np.append(observablesErr_out, spacing_errs)
        filename_suffix+='-dP'
        observables.remove('period_spacing')

    if 'rope_length' in observables:
        RL, error = PdP_pattern_rope_length(periods, P_error=periodsErr)
        observables_out = np.append(observables_out, RL)      # Include pattern rope length as observable
        observablesErr_out = np.append(observablesErr_out, error)
        filename_suffix+='-rl'
        observables.remove('rope_length')

    # Add all other observables in the list from the dataFrame
    for observable in observables:
        Obs = np.asarray(Obs_dFrame[observable]) # since these columns have less entries, the dataFrame has NaNs in the empty rows
        Obs = Obs[~np.isnan(Obs)]                # remove all entries that are NaN in the numpy array
        ObsErr = np.asarray(Obs_dFrame[f'{observable}_err'])
        ObsErr = ObsErr[~np.isnan(ObsErr)]

        observables_out = np.append(observables_out, Obs)
        observablesErr_out = np.append(observablesErr_out, ObsErr)
        filename_suffix+=f'-{observable}'

    return observables_out, observablesErr_out, filename_suffix

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def mle_chi2(YObs, ObsErr, YTheo):
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
def mle_mahalanobis(YObs, ObsErr, YTheo):
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
    check_matrix(V)     # check if positive definite
    # Calculate Mahalanobis distances
    MD = np.zeros(q)
    Vinv = np.linalg.inv(V)
    for i in range(q):
        diff = (YTheoMat[:,i]-YObsMat)
        MD[i] = np.matmul(np.matmul(diff.T,Vinv),diff)[0][0]

    return MD

################################################################################
def check_matrix(V, plot=False):
    """
    Check the if the the eigenvalues of the Variance matrix are all positive,
    since this means the matrix is positive definite.
    ------- Parameters -------
    V: 2D np array
        Variance matrix
    plot: boolean
        Flag to show a plot of the variance matrix
    """
    if np.all(np.linalg.eigvals(V) > 0)==False: # If all eigencalues are >0, it is positive definite
        sys.exit(logger.error('V matrix is possibly not positive definite (since eigenvalues are not all > 0)'))

    if plot is True:
        im = plt.imshow(V, aspect='auto')
        plt.colorbar(im)
        plt.show()

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def mle_rope_length(obs_P, obs_P_error, theo_P):
    """
    Calculate rope length values for the given theoretial patterns. So for each pattern
    the total sum of distances between the theoretical points and their observed error boxes in the P-dP plane.
    ------- Parameters -------
    YObs, ObsErr: numpy array of floats
        Observed values and their errors (period or frequency in units of day or 1/day )
    YTheo: numpy array of arrays of floats
        Array of all theoretical patterns to calculate the chi squared value for.

    ------- Returns -------
    rope length: numpy array of floats
        rope length values for the given theoretical values
    """
    rope_length = []
    obs_dP, obs_dP_error = ffg.generate_obs_series(obs_P, obs_P_error)  # observed period spacing in seconds
    obs_dP = np.asarray(obs_dP)#/86400     TODO     # switch back from seconds to days (so both P and dP are in days)

    for i in range(len(theo_P)):
        theo_dP = ffg.generate_thry_series(theo_P[i])   # theoretical period spacing in seconds
        theo_dP = np.asarray(theo_dP)#/86400  TODO  # switch back from seconds to days (so both P and dP are in days)

        total_length = 0
        for j in range(len(obs_dP)):
            total_length += distance_point_errorBox( obs_P[j], obs_dP[j], obs_P_error[j], obs_dP_error[j], theo_P[i][j], theo_dP[j])

        rope_length.append(total_length)
    rope_length = np.asarray(rope_length)

    return rope_length
################################################################################
def distance_point_errorBox(obs_P, obs_dP, obs_P_error, obs_dP_error, theo_P, theo_dP):
    """
    Calculate the shortest distance from a theoretical point to the observed error box. Distance 0 if it is in the error box.
    ------- Parameters -------
    obs_P, obs_dP, theo_P, theo_dP: float
        Observed (obs) and theoretical (theo) periods (P) and period spacing (dP)
    obs_P_error, obs_dP_error: float
        Error on the observed period and period spacing

    ------- Returns -------
    float
        The distance between the theoretical point and the observed error box.
    """
    # calculate edges fo the error box
    box_x_min = obs_P - obs_P_error
    box_x_max = obs_P + obs_P_error
    box_y_min = obs_dP - obs_dP_error
    box_y_max = obs_dP + obs_dP_error

    # Look at the tuple being provided to the max function: (min-p, 0, p-max). Let's designate this tuple (a,b,c).
    # If p is left of min, then we have p < min < max, which means the tuple will evaluate to (+,0,-), and so the max function will correctly return a = min - p.
    # If p is between min and max, then we have min < p < max, which means the tuple will evaluate to (-,0,-). So again, the max function will correctly return b = 0.
    # Lastly, if p is to the right of max, then we have, min < max < p, and the tuple evaluates to (-,0,+). Once again, max correctly returns c = p - max.
    dx = np.max([box_x_min - theo_P,  0, theo_P - box_x_max])
    dy = np.max([box_y_min - theo_dP, 0, theo_dP - box_y_max])
    return np.sqrt(dx**2+dy**2)

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

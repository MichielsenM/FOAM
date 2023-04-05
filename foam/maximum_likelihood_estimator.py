""" Functions to perform different kinds of maximum likelihood estimation for the models in a grid, and make correlation plots.
Note: The file with observations needs to hold temperature as Teff, although the analysis is done using the logTeff values."""
# from foam import maximum_likelihood_estimator as mle
import numpy as np
import pandas as pd
import sys, logging, os
import matplotlib.pyplot as plt
from pathlib import Path
from foam import support_functions as sf
from foam import functions_for_gyre as ffg
from foam import functions_for_mesa as ffm

logger = logging.getLogger('logger.mle_estimator')  # Make a child logger of "logger" made in the top level script

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def calculate_likelihood(Obs_path, Theo_file, observables=None, merit_function=None, star_name=None, fixed_params=None):
    """
    Perform a maximum likelihood estimation using the provided type of merit function on the list of  observables.
    Writes a data file with the values of the merit funtion and input parameters of each model.
    Can also select and continue the analysis of nested grids through the keyword 'fixed_params'.
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
    fixed_params: dictionary
        Only select and analyse the part of the theoretical grid with the specified parameter values.
        The keys specify for which parameters only the specified value should be selected.
    """
    # Read in the observed data and make an array of the observed obervables
    Obs_dFrame = pd.read_table(Obs_path, delim_whitespace=True, header=0)
    Obs, ObsErr, file_suffix_observables = create_obs_observables_array(Obs_dFrame, observables)

    Path_theo   = Path(Theo_file)
    #suffix for filename to indicate the merit funtion used
    suffix = {'chi2'       : 'CS',
              'mahalanobis': 'MD'}

    # set the name of the output file and make it's directory if needed
    head, tail = sf.split_line(Path_theo.stem, star_name)
    DataOutDir = Path(f'{os.getcwd()}/{str(Path_theo.parent).split("/")[-1]}')
    Path(DataOutDir).mkdir(parents=True, exist_ok=True)
    DataOut = f'{DataOutDir}/{star_name}{tail}_{suffix[merit_function]}_{file_suffix_observables}.hdf'

    # Theoretical grid data
    Theo_dFrame = sf.get_subgrid_dataframe(Theo_file,fixed_params)
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
    if len(neg_value)>0:
        logger.warning(f"""{len(neg_value)} models were discarded due to mismatches during the selection of theoretical frequencies.
                        This is likely due to the frequency range used for the theoretical calculations being too narrow.""")
    Theo_observables = np.asarray(newTheo)
    Thetas           = np.asarray(newThetas)

    # Dictionary containing different merit functions
    switcher={ 'chi2': merit_chi2,
                'mahalanobis' : merit_mahalanobis}

    # get the desired function from the dictionary. Returns the lambda function if option is not in the dictionary.
    selected_merit_function = switcher.get(merit_function, lambda x, y, z: sys.exit(logger.error('invalid type of maximum likelihood estimator')))
    merit_values = selected_merit_function(Obs, ObsErr, Theo_observables, fig_title=f'{star_name}{tail}_{suffix[merit_function]}_{file_suffix_observables}', star_name=star_name)

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
    df.to_hdf(f'{DataOut}', 'merit_values', format='table', mode='w')
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
            spacing, _ = ffg.generate_spacing_series(periods_part)
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
            spacing, spacing_errs = ffg.generate_spacing_series(periods, periodsErr)
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
def merit_chi2(YObs, ObsErr, YTheo, fig_title=None, star_name=None):
    """
    Calculate chi squared values for the given theoretical patterns
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
def merit_mahalanobis(YObs, ObsErr, YTheo, fig_title=None, star_name=None):
    """
    Calculate mahalanobis distance (MD) values for the given theoretical patterns.
    ------- Parameters -------
    YObs, ObsErr: numpy array of floats
        Observed values and their errors (period or frequency)
    YTheo: numpy array of arrays of floats
        Array of all theoretical patterns to calculate the MD value for.

    ------- Returns -------
    MD: numpy array of floats
        Mahalanobis distances for the given theoretical patterns.
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
    check_matrix(V, fig_title=fig_title, star_name=star_name)     # check if positive definite and make figure
    # Calculate Mahalanobis distances
    MD = np.zeros(q)
    Vinv = np.linalg.inv(V)
    for i in range(q):
        diff = (YTheoMat[:,i]-YObsMat)
        MD[i] = np.matmul(np.matmul(diff.T,Vinv),diff)[0][0]

    return MD

################################################################################
def check_matrix(V, plot=True, fig_title='Vmatrix', star_name=None):
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
    file_Path = Path(f'{os.getcwd()}/V_matrix/{star_name}_determinant_conditionNr.tsv')
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
        plt.clf()
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
        dP, dP_error = ffg.generate_spacing_series(P, P_error)
    else:
        dP, _ = ffg.generate_spacing_series(P)

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

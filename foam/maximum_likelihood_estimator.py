""" Functions to perform different kinds of maximum likelihood estimation for the models in a grid, and make correlation plots.
Note: The file with observations needs to hold temperature as Teff, although the analysis is done using the logTeff values."""
import numpy as np
import pandas as pd
import sys, logging, os
import matplotlib.pyplot as plt
from pathlib import Path
from foam import support_functions as sf
from foam import build_optimised_pattern as bop

logger = logging.getLogger('logger.mle_estimator')  # Make a child logger of "logger" made in the top level script

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def calculate_likelihood(Theo_file, observables=None, merit_function=None, Obs_path=None, star_name=None, fixed_params=None, grid_parameters=None):
    """
    Perform a maximum likelihood estimation using the provided type of merit function on the list of observables.
    Writes a data file with the values of the merit funtion and input parameters of each model.
    Can select and continue the analysis of nested grids through the keyword 'fixed_params'.
    ------- Parameters -------
    Obs_path: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    Theo_file: string
        Path to the hdf5 file with the theoretical model input parameters (first set of columns), frequency or period values (last set of columns),
        and possibly extra columns with additional observables (these columns should be in between the input parameters and frequency columns).
    observables: list of strings
        Which observables are included in the merit function.
        Must contain either 'f' (frequency), 'P' (period), or 'dP' (period-spacing) which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.
    merit_function: string
        The type of merit function to use. Currently supports "CS and "MD" ("chi-squared" and "mahalanobis distance").
    star_name: string
        Name of the star, used in file naming.
    fixed_params: dictionary
        Only select and analyse the part of the theoretical grid with the specified parameter values.
        The keys specify for which parameters only the specified value should be selected.
    grid_parameters: list of string
        List of the parameters in the theoretical grid.
    """
    if 'f' in observables:
        observed_quantity = 'frequency'
    elif 'P' or 'dP' in observables:
        observed_quantity = 'period'
    # Read in the observed data and make an array of the observed obervables
    Obs_dFrame = pd.read_table(Obs_path, delim_whitespace=True, header=0, index_col='index')
    Obs, ObsErr, file_suffix_observables = create_obs_observables_array(Obs_dFrame, observables)

    # set the name of the output file
    _, tail = sf.split_line(Path(Theo_file).stem, star_name)
    filename = f'{star_name}{tail}_{merit_function}_{file_suffix_observables}'

    # Theoretical grid data
    Theo_dFrame = sf.get_subgrid_dataframe(Theo_file,fixed_params)

    missing_absolute = np.where(Theo_dFrame.columns.to_series().str.contains(f'{observed_quantity}_missing'))[0]    # get the interruptions in the pattern, absolute index in dataframe
    missing_relative = np.where(Theo_dFrame.filter(like=f'{observed_quantity}').columns.to_series().str.contains(f'{observed_quantity}_missing'))[0]  # get the interruptions in the pattern, index relative within pulsations
    Theo_dFrame = Theo_dFrame.drop(columns=Theo_dFrame.columns[missing_absolute])   # Remove columns of missing frequencies
    missing_indices=[ missing_relative[i]-i for i in range(len(missing_relative)) ] # Adjust indices for removed lines of missing frequencies
    
    Thetas    = np.asarray(Theo_dFrame.filter(['rot']+['rot_err']+grid_parameters ))
    Theo_puls = np.asarray(Theo_dFrame.filter(like=f'{observed_quantity}'))

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
    newTheo = np.asarray(newTheo)
    newThetas = np.asarray(newThetas)

    # Dictionary containing different merit functions
    switcher={'CS': merit_chi2,
              'MD' : merit_mahalanobis}

    # get the desired function from the dictionary. Returns the lambda function if option is not in the dictionary.
    selected_merit_function = switcher.get(merit_function, lambda x, y, z: sys.exit(logger.error('invalid type of maximum likelihood estimator')))
    merit_values = selected_merit_function(Obs, ObsErr, newTheo, fig_title=f'{filename}', star_name=star_name)

    # Combine values and save the results
    CombData = np.concatenate((np.matrix(merit_values).T,newThetas),axis=1)  # add an additional column for MLE 'meritValues'
    df = pd.DataFrame(data=CombData, columns=['meritValue']+['rot']+['rot_err']+grid_parameters) # put the data in a pandas DataFrame
    df = pd.merge(df, Theo_dFrame.drop(list(Theo_dFrame.filter(regex=f'{observed_quantity}').columns), axis=1), how='inner', on=['rot', 'rot_err']+grid_parameters)
    df.to_hdf(f'{os.getcwd()}/meritvalues/{filename}.hdf', 'merit_values', format='table', mode='w')
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
        Must contain either 'f' (frequency), 'P' (period), or 'dP' (period-spacing) which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.
    missing_indices: list of int
        Contains the indices of the missing pulsations so that the period sapcing pattern can be split around them.
    ------- Returns -------
    observables_out: numpy array of floats
        The values of the specified observables for the model.
    """
    observables = list(observables_in)  # Make a copy to leave the array handed to this function unaltered.

    if 'P' in observables:   # Cast to list before using asarray to prevent memory leak
        observables_out = np.asarray(list(Theo_dFrame.filter(like='period').loc[index]))  # add the periods to the output list
        observables.remove('P')

    elif 'f' in observables: # Cast to list before using asarray to prevent memory leak
        observables_out = np.asarray(list(Theo_dFrame.filter(like='frequency').loc[index]))  # add the frequencies to the output list
        observables.remove('f')

    elif 'dP' in observables:
        periods = np.asarray(Theo_dFrame.filter(like='period').loc[index])
        observables_out = []
        for periods_part in np.split(periods,missing_indices):
            spacing, _ = bop.generate_spacing_series(periods_part)
            observables_out = np.append(observables_out, np.asarray(spacing)/86400) # switch back from seconds to days (so both P and dP are in days)
        observables.remove('dP')

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
        Which observables are included in the returned array.
        Must contain either 'f' (frequency), 'P' (period), or 'dP' (period-spacing) which will be computed for the period pattern.
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
    filename_suffix = ''

    if 'P' in observables:
        observables_out = np.asarray(Obs_dFrame['period'])
        observablesErr_out = np.asarray(Obs_dFrame['period_err'])
        filename_suffix = 'P'
        observables.remove('P')

    elif 'f' in observables:
        observables_out = np.asarray(Obs_dFrame['frequency'])
        observablesErr_out = np.asarray(Obs_dFrame['frequency_err'])
        filename_suffix = 'f'
        observables.remove('f')

    elif 'dP' in observables:
        observables_out = []
        observablesErr_out = []
        period = np.asarray(Obs_dFrame['period'])
        periodErr = np.asarray(Obs_dFrame['period_err'])
        periods_parts = np.split(period,missing_indices)
        periodsErr_parts = np.split(periodErr,missing_indices)
        for periods, periodsErr in zip(periods_parts, periodsErr_parts):
            spacing, spacing_errs = bop.generate_spacing_series(periods, periodsErr)
            observables_out = np.append(observables_out, np.asarray(spacing)/86400) # switch back from seconds to days (so both P and dP are in days)
            observablesErr_out = np.append(observablesErr_out, np.asarray(spacing_errs)/86400)

        filename_suffix = 'dP'
        observables.remove('dP')

    if len(observables) > 0:
        filename_suffix+='+extra'
    # Add all other observables in the list from the dataFrame
    for observable in observables:
        logTeff=False
        if observable == 'logTeff': 
            observable = 'Teff' # To read it as Teff from the observations datafile
            logTeff = True
        Obs = np.asarray(Obs_dFrame[observable]) # since these columns have less entries, the dataFrame has NaNs in the empty rows
        Obs = Obs[~np.isnan(Obs)]                # remove all entries that are NaN in the numpy array
        ObsErr = np.asarray(Obs_dFrame[f'{observable}_err'])
        ObsErr = ObsErr[~np.isnan(ObsErr)]

        if logTeff: # Convert observed Teff and error to log values
            logT = np.log10(Obs)
            logT_err = ObsErr / (Obs * np.log(10))
            observables_out = np.append(observables_out, logT)
            observablesErr_out = np.append(observablesErr_out, logT_err)
        else:
            observables_out = np.append(observables_out, Obs)
            observablesErr_out = np.append(observablesErr_out, ObsErr)

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
    fig_title, star_name: None
        Should not be used in this function, but is to make it analogous to merit_mahalanobis()
        and enable the use of the lambda function.
    ------- Returns -------
    chi2: numpy array of floats
        chi squared values for the given theoretical values
    """
    chi2 = np.array([  np.sum(((one_YTheo-YObs)/ObsErr)**2) for one_YTheo in YTheo ])
    return chi2

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def merit_mahalanobis(YObs, ObsErr, YTheo, generate_output=True, fig_title=None, star_name=None):
    """
    Calculate mahalanobis distance (MD) values for the given theoretical patterns.
    ------- Parameters -------
    YObs, ObsErr: numpy array of floats
        Observed values and their errors (period or frequency)
    YTheo: numpy array of arrays of floats
        Array of all theoretical patterns to calculate the MD value for.
    generate_output: boolean
        Flag to write output and plot the variance-covariance matrix
    fig_title: string
        The name of the figure to be created.
    star_name: string
        The name of the analysed star, for file naming purposes.
    ------- Returns -------
    MD: numpy array of floats
        Mahalanobis distances for the given theoretical patterns.
    """
    # Convert to matrix format (np.matrix is not recommended, use array and nexaxis instead)
    YObsMat = np.array(YObs)[np.newaxis].T
    YTheoMat = np.array(YTheo)[np.newaxis].T

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
    check_matrix(V, generate_output=generate_output, fig_title=fig_title, star_name=star_name)     # check if positive definite and make figure
    # Calculate Mahalanobis distances
    MD = np.zeros(q)
    Vinv = np.linalg.inv(V)
    for i in range(q):
        diff = (YTheoMat[:,i]-YObsMat)
        MD[i] = np.matmul(np.matmul(diff.T,Vinv),diff)[0][0]

    return MD

################################################################################
def check_matrix(V, generate_output=True, fig_title='Vmatrix', star_name=None):
    """
    Check the if the the eigenvalues of the Variance-covariance matrix are all positive,
    since this means the matrix is positive definite. Compute its determinant and condition number,
    and write them to a tsv file. Create and save a figure of the variance-covariance matrix.
    ------- Parameters -------
    V: 2D np array
        Variance-covariance matrix
    output: boolean
        Flag to write output and plot the variance-covariance matrix
    fig_title: string
        The name of the figure to be created.
    star_name: string
        The name of the analysed star, for file naming purposes.
    """
    if np.all(np.linalg.eigvals(V) > 0)==False: # If all eigenvalues are >0, it is positive definite
        logger.error('V matrix is possibly not positive definite (since eigenvalues are not all > 0)')
        sys.exit(1)
    logger.debug(f'max(V) = {np.max(V)}')
    kk=10 # multiply the matrix by the exponent of this, otherwise the determinant can be too small for the numerics

    if generate_output is True:
        file_Path = Path(f'{os.getcwd()}/V_matrix/{star_name}_determinant_conditionNr.tsv')
        file_Path.parent.mkdir(parents=True, exist_ok=True)
        if not file_Path.is_file():
            with file_Path.open("w") as file:
                file.write(f'method \t ln(det(V)) \t condition_number \n')
        with file_Path.open("a") as file:
            file.write(f'{fig_title} \t {np.log(np.linalg.det(np.exp(kk)*V))-kk*V.shape[0]:.2f} \t {np.linalg.cond(V):.2f} \n ')


        im = plt.imshow(V*10**4, aspect='auto', cmap='Reds') # Do *10^4 to get rid of small values, and put this in the colorbar label
        plt.ylabel(rf'Obs {V.shape[0]}  $\leftarrow$ Obs 1', size=14)
        plt.xlabel(rf'Obs 1 $\rightarrow$ Obs {V.shape[0]} ', size=14)

        cbar = plt.colorbar(im)
        cbar.ax.set_ylabel(r'[d$^{2} 10^{-4}$]', rotation=90, labelpad=15, size=14)
        cbar.ax.tick_params(labelsize=14)
        plt.title('Variance covariance matrix')
        plt.tick_params(labelsize=14)
        plt.tight_layout()
        plt.savefig(f'{os.getcwd()}/V_matrix/{fig_title}.png')
        plt.clf()
        plt.close('all')

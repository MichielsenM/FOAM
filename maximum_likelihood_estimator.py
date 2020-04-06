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

def plot_correlations(MLE_values_file, MESA_grid_dir = '/lhome/mathiasm/MESA_grid_ECP-DE', fig_title=None, label_size=18, fig_outputDir='figures/',
                      percentile_to_show=0.1, logTeff_obs=[np.log10(11650), [0.00790, 0.00776]], logg_obs=[3.97, 0.08]):
    """
    Make a plot of all variables vs each other variable (and also a Kiel/HRD diagram), showing the MLE values as colorscale.
    The subplots on the diagonal show the distribution of that variable.
    The list of variables is retrieved from columns of the MLE_values_file,
    where the first column is 'distance', which are the MLE values.
    The resulting figure is saved afterwards in the specified location.
    ------- Parameters -------
    MLE_values_file, MESA_grid_dir: string
        path to the file with the MLE values and to the MESA profiles directory
    fig_title: string
        Title of the figure
    label_size: int
        size of the axis labels
    fig_outputDir: string
        output directory for the figures
    percentile_to_show: float
        percentile of models to show in the plots
    logTeff_obs: list (length 2)
        log of the observed effective temperature (first entry) and its error (second entry, again a list for non-symmetric errors)
    logg_obs: list (length 2) of float
        log of the observed surface gravity (first entry) and its error (second entry)
    """
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
    ax_hrd.set_xlabel('log Teff')
    ax_hrd.set_ylabel('log g')

    Teff = []
    # logL = []
    logg = []
    color  = []
    if 'ECP_' in MLE_values_file:
        CBM = 'extended_convective_penetration'
    elif 'DO_' in MLE_values_file:
        CBM = 'diffusive_overshoot'

    for i in range(0, df.shape[0]):
        prof =  f'{MESA_grid_dir}/{CBM}/MESA_out/Zini{df.iloc[i]["Z"]}/profiles/'    \
                +f'Z{df.iloc[i]["Z"]}_M{df.iloc[i]["M"]:.2f}_logD{df.iloc[i]["logD"]:.2f}_'    \
                +f'aov{df.iloc[i]["aov"]:.3f}_fov{df.iloc[i]["fov"]:.3f}_Xc{df.iloc[i]["Xc"]:.2f}.h5_prof'

        prof_data = mypy.read_hdf5(prof)
        Teff.append(np.log10(float(prof_data['Teff'])))
        logg.append(prof_data['log_g'][0])
        color.append(np.log10(df.iloc[i,0]))
        # logL.append(np.log10(float(prof_data['photosphere_L'])))

    im = ax_hrd.scatter(Teff, logg, c=color, cmap='hot')
    ax_hrd.errorbar(logTeff_obs[0], logg_obs[0], xerr=np.array([logTeff_obs[1]]).T, yerr=logg_obs[1])    # Spectroscopic error bar

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
    cbar.set_label('log(Mahalanobis distance)', rotation=90)
    plt.subplots_adjust(left=0.11, right=0.87, bottom=0.1, top=0.95)

    fig.suptitle(fig_title)
    Path(fig_outputDir).mkdir(parents=True, exist_ok=True)
    fig.savefig(f'{fig_outputDir}{fig_title}.png', dpi=300)
    plt.close('all')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def estimate_max_likelihood(Obs_path, Theo_file, which_observable='', estimator_type ='', dP_in_Vmatrix=False):
    """
    Perform a maximum likelihood estimation using the provided specifications (type of estimator, observables ...)
    Writes a data file with the MLE values and input parameters of each model.
    ------- Parameters -------
    Obs_path: string
        Path to the tsv file with observations, which has the following format: row 1 header,
        row 2 and 3 observed frequencies and their errors, row 4 and 5 observed periods and their errors.
    Theo_file: string
        Path to the tsv file with the theoretical frequency or period values and model input parameters.
    which_observable: string
        Which observables are used in the analysis, options are 'frequencies' or 'periods'.
    estimator_type: string
        The type of maximum likelihood estimation to use. Currently supports "chi2", "mahalanobis" and "rope_length"
    dP_in_Vmatrix: boolean
        If False, make the V matrix (for mahalanobis distances) with only frequency or period values
        if True, (observable should be periods) include both periods and period spacing values in the V matrix
    """
    if dP_in_Vmatrix is True and (which_observable != 'periods' or estimator_type != 'mahalanobis'):
        logger.warning('Skipped this function call, period spacing in variance matrix must use periods and mahalanobis distances.'\
                        +f' However which_observable={which_observable} and estimator_type={estimator_type}')
        return

    # Read in the desired observables
    ObsVal = pd.read_table(Obs_path, delim_whitespace=True, header=0)
    if which_observable == 'frequencies':
        Obs = np.asarray(ObsVal.iloc[0])
        ObsErr = np.asarray(ObsVal.iloc[1]*4)   # times 4 to take correlation structure of the data into account (see Moravveji et al. 2016)
    elif which_observable == 'periods':
        Obs = np.asarray(ObsVal.iloc[2])
        ObsErr = np.asarray(ObsVal.iloc[3]*4)   # times 4 to take correlation structure of the data into account (see Moravveji et al. 2016)
    else:
        sys.exit(logger.error(f'Observable \"{which_observable}\" not recognised, should be \"periods\" or \"frequencies\"' ))

    # Grid data
    TheoFile = pd.read_table(Theo_file, delim_whitespace=True, header=0)
    Thetas =  np.asarray(TheoFile.loc[:,:'Xc']) # varied parameters in the grid (e.g. Mini, Xini, Xc etc.)
    Theo =  np.asarray(TheoFile.loc[:,'f1':])   # corresponding theoretical values to the observed ones
    Path_theo = Path(Theo_file)

    # Make new list of theoretical models without entries with value -1
    newTheo = []
    newThetas = []
    for i in range(len(Theo)):
        if (Theo[i][0] !=-1) and (Theo[i][-1] !=-1):     # ignore models where one of the freqs is -1
            if dP_in_Vmatrix is True:
                spacing = ffg.generate_thry_series(Theo[i])
                spacing = np.asarray(spacing)/86400 # switch back from seconds to days (so both P and dP are in days)
                aux = np.append(Theo[i], spacing)   # Include both P and dP as observables
            else:
                aux = Theo[i]

            newTheo.append(aux)
            newThetas.append(Thetas[i])
    neg_value = np.unique(np.where(Theo==-1)[0])

    logger.debug(f'ignored -1 freqs : {len(neg_value)}')
    logger.debug(f'original #models : {len(Theo)}')
    logger.debug(f'remaining #models: {len(newTheo)}')

    Theo = np.asarray(newTheo)
    Thetas = np.asarray(newThetas)

    #suffix for filename to indicate the MLE method used
    suffix = {'chi2'       : 'CS',
              'mahalanobis': 'MD',
              'rope_length': 'RL'}

    # Add period spacing to the list of observables
    if dP_in_Vmatrix is True:
        observed_spacings,observed_spacings_errors = ffg.generate_obs_series(Obs, ObsErr)

        observed_spacings = np.asarray(observed_spacings)/86400 # switch back from seconds to days (so both P and dP are in days)
        observed_spacings_errors = np.asarray(observed_spacings_errors)/86400   # switch back from seconds to days (so both P and dP are in days)

        Obs = np.append(Obs, observed_spacings) # Include both P and dP as observables
        ObsErr = np.append(ObsErr, observed_spacings_errors)

        head, tail = mypy.split_line(Path_theo.stem, '_')
        DataOut = f'{Path_theo.parent}/PdP_{tail}_{suffix[estimator_type]}.dat'
    else:
        DataOut = f'{Path_theo.parent}/{Path_theo.stem}_{suffix[estimator_type]}.dat'

    # Dictionary containing different functions for MLE
    switcher={ 'chi2': mle_chi2,
                'mahalanobis' : mle_mahalanobis,
                'rope_length': mle_rope_length }

    # get the desired function from the dictionary. Returns the lambda function if option is not in the dictionary.
    mle_function = switcher.get(estimator_type, lambda x, y, z: sys.exit(logger.error('invalid type of maximum likelihood estimator')))
    mle_values = mle_function(Obs, ObsErr, Theo)

    # Print smallest and highest values
    idx2 = np.argsort(mle_values)
    Parameters = TheoFile.loc[:,:'Xc'].columns  # Parameter names
    logger.info(f'Smallest {estimator_type} distance: {mle_values[np.argsort(mle_values)][0]}')
    logger.info(f'Highest {estimator_type} distance : {mle_values[np.argsort(mle_values)][-1]}')
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

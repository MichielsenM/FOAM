"""Helpful functions for GYRE input and output, (e.g. extracting all frequencies in a grid to 1 file,
   constructing theoretical pulsation patterns, calculate GYRE scanning range to find desired radial orders...)"""
# from foam import functions_for_gyre as ffg
import numpy as np
import pandas as pd
from pathlib import Path
import glob, sys, pkgutil
import logging
import multiprocessing
from functools import partial
from io import StringIO
from foam import support_functions as sf
from foam import functions_for_mesa as ffm

logger = logging.getLogger('logger.ffg')

################################################################################
def generate_spacing_series(periods, errors=None):
    """
    Generate the period spacing series (delta P = p_(n+1) - p_n )
    ------- Parameters -------
    periods, errors (optional): list of floats
        Periods and their errors in units of days
    ------- Returns -------
    observed_spacings, observed_spacings_errors: list of floats
        period spacing series (delta P values) and its errors (if supplied) in units of seconds
    """
    spacings = []
    if errors is None:
        spacings_errors = None
        for n,period_n in enumerate(periods[:-1]):
            spacings.append( abs(period_n - periods[n+1])*86400. )
    else:
        spacings_errors = []
        for n,period_n in enumerate(periods[:-1]):
            spacings.append( abs(period_n - periods[n+1])*86400. )
            spacings_errors.append(np.sqrt( errors[n]**2 + errors[n+1]**2 )*86400.)

    return spacings, spacings_errors

################################################################################
def extract_frequency_grid(gyre_files, output_file='pulsationGrid.tsv', parameters=['rot', 'Z', 'M', 'logD', 'aov', 'fov', 'Xc']):
    """
    Extract frequencies from each globbed GYRE file and write them to 1 large file.
    ------- Parameters -------
    gyre_files: string
        String to glob to find all the relevant GYRE summary files.
    output_file: string
        Name (can include a path) for the file containing all the pulsation frequencies of the grid.
    parameters: list of strings
        List of parameters varied in the computed grid, so these are taken from the
        name of the summary files, and included in the 1 file containing all the info of the whole grid.
    """
    # make a copy of the list, so parameters is not extended with all the orders before passing it on to 'all_freqs_from_summary'
    header_parameters = list(parameters)

    df = pd.DataFrame(columns=header_parameters)  # dataframe for all the pulations
    MP_list = multiprocessing.Manager().list()    # make empty MultiProcessing listProxy

    # Glob all the files, then iteratively send them to a pool of processors
    summary_files = glob.iglob(gyre_files)
    p = multiprocessing.Pool()
    extract_func = partial(all_freqs_from_summary, parameters=parameters)
    dictionaries = p.imap(extract_func, summary_files)
    for new_row in dictionaries:
        MP_list.append(new_row)   # Fill the listProxy with dictionaries for each read file

    df = df.append(MP_list[:], ignore_index=True) # Combine the dictionaries into one dataframe
    # Sort the columns with frequencies by their radial order
    column_list = list(df.columns[:len(parameters)])
    column_list.extend(sorted(df.columns[len(parameters):]))
    df = df.reindex(column_list, axis=1)

    # Generate the directory for the output file and write the file afterwards
    Path(Path(output_file).parent).mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, sep='\t',index=False) # write the dataframe to a tsv file
    p.close()
################################################################################
def all_freqs_from_summary(GYRE_summary_file, parameters):
    """
    Extract model parameters and pulsation frequencies from a GYRE summary file
    ------- Parameters -------
    GYRE_summary_file: string
        path to the GYRE summary file
    parameters: list of strings
        List of input parameters varied in the computed grid, so these are read from the filename and included in returned line.

    ------- Returns -------
    param_dict: dictionary
        Dictionary containing all the model parameters and pulsation frequencies of the GYRE summary file.
    """

    data = sf.read_hdf5(GYRE_summary_file)
    param_dict = sf.get_param_from_filename(GYRE_summary_file, parameters)

    for j in range(len(data['freq'])-1, -1, -1):    # Arrange increasing in radial order
        n_pg = data["n_pg"][j]
        if abs(n_pg) < 10:
            n_pg = f'{sf.sign(n_pg)}0{abs(n_pg)}'
        param_dict.update({f'n_pg{n_pg}':data['freq'][j][0]})
        # param_dict.update({f'n_pg{data["n_pg"][j]}':data['freq'][j][0]})

    return param_dict

################################################################################
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
def ledoux_splitting(frequencies, betas, Mstar, Rstar, omega=0, m=1):
    """
    Calculate rotationally shifted frequencies in a perturbative way. (See Aerts et al. (2010), eq 3.357)
    ------- Parameters -------
    frequencies, betas: numpy array of floats
        frequency values (c/d) and beta values (eq 3.357, Aerts et al. (2010))
    Mstar, Rstar, omega: float
        stellar mass (g), radius (cm) and rotation frequency in units of critical velocity (omega_crit^-1)
    m: int
        azimuthal order

    ------- Returns -------
    shifted_freqs: numpy array of floats
        pulsation frequency values, shifted by the Ledoux splitting
    """

    G = 6.67428E-8 # gravitational constant (g^-1 cm^3 s^-2)
    omega_crit = (1/(2*np.pi))*(8*G*Mstar/(27*Rstar**3))**0.5 # Roche critical rotation frequency (s^-1)
    omega_cycday = omega*omega_crit*86400 # rotation frequency in units of c/d

    shifted_freqs = frequencies-(m*omega_cycday*(1-betas)) # shifted frequencies in units of c/d
    return shifted_freqs

################################################################################
def calc_scanning_range(gyre_file_path, npg_min=-50, npg_max=-1, l=1, m=1, omega_rot=0.0, unit_rot = 'CYC_PER_DAY', rotation_frame='INERTIAL'):
    """
    Calculate the frequency range for the sought radial orders of the g modes.
    ------- Parameters -------
    gyre_file_path: string
        absolute path to the gyre file that needs to be scanned
    n_min, n_max: integer
        lower and upper values of the required range in radial order
    l, m: integer
        degree (l) and azimuthal order (m) of the modes
    omega_rot: float
        rotation frequency of the model
    unit_rot: string
        unit of the rotation frequency, can be CYC_PER_DAY or CRITICAL (roche critical)
    rotation_frame: string
        rotational frame of reference for the pulsation freqencies

    ------- Returns -------
    f_min, f_max: float
        lower and upper bound of frequency range that needs to be scanned in oder
        to retrieve the required range of radial orders
    """
    directory, gyre_file = sf.split_line(gyre_file_path, 'gyre/') # get directory name and GYRE filename
    Xc_file = float(sf.substring(gyre_file, 'Xc', '.GYRE'))       # get Xc
    MESA_hist_name, tail = sf.split_line(gyre_file, '_Xc')        # Get the MESA history name form the GYRE filename
    hist_file = glob.glob(f'{directory}history/{MESA_hist_name}.*hist')[0]   # selects MESA history file (.hist or .h5_hist) corresponding to the GYRE file

    header, data  = ffm.read_mesa_file(hist_file)
    Xc_values = np.asarray(data['center_h1'])
    P0_values = np.asarray(data['Asymptotic_dP'])

    # Obtain the asymptotic period spacing value/buoyancy radius at the Xc value closest to that of the gyre file
    diff = abs(Xc_file - Xc_values)
    xc_index = np.where(diff == np.min(diff))[0]
    P0 = P0_values[xc_index][0]/86400 # asymptotic period spacing value/buoyancy radius, /86400 to go from sec to day

    # Calculate the scanning range a bit broader than the purely asymptotic values, just to be safe.
    n_max_used = abs(npg_min-3)
    n_min_used = abs(min(-1, npg_max+3))

    if omega_rot==0:
        # If no rotation, use asymptotic values
        f_min = np.sqrt(l*(l+1)) / (n_max_used*P0)
        f_max = np.sqrt(l*(l+1)) / (n_min_used*P0)
    else:
        if unit_rot == 'CRITICAL': # Roche critical
            model_mass   = ffm.convert_units('mass',   np.asarray(data['star_mass'])[xc_index], convertto='cgs')
            model_radius = ffm.convert_units('radius', 10**np.asarray(data['log_R'])[xc_index], convertto='cgs')
            G = ffm.convert_units('cgrav', 1)
            Roche_rate = (1/(2*np.pi))*np.sqrt((8*G*model_mass)/(27*model_radius**3)) # Roche crit rotation rate in cycles per second
            Roche_rate = Roche_rate * 86400 # Roche crit rotation rate in cycles per day
            omega_rot = omega_rot * Roche_rate # Multiply by fraction of the crit rate, to get final omega_rot in cycles per day

        # Make a pandas dataframe containing an interpolation table for lambda (eigenvalues of LTE - TAR)
        data = pkgutil.get_data(__name__, 'lambda.csv')
        data_io = StringIO(data.decode(sys.stdout.encoding))
        df = pd.read_csv(data_io, sep=",")

        # will add extra functionality to calculate the bounds explicitly, making use of GYRE

        # Select nu (spin parameter) and lambda column when values in l and m column correspond to requested values
        ###### SHOULD BE CHANGED TO PARAMETER 'K' ---> needs adjustment in lambda.csv - JVB.
        NuLambda = df.loc[(df['l'] == l) & (df['m'] == m)][['nu', 'Lambda']]

        # Generate numpy array from pandas dataframe series
        nu = NuLambda['nu'].to_numpy()
        Lambda = NuLambda['Lambda'].to_numpy()

        # Generate difference between pulsation frequency and asymptotic value (in co-rotating frame) in units of c/d
        diff_max = nu/(2.*omega_rot) - P0*n_max_used/np.sqrt(Lambda)
        diff_min = nu/(2.*omega_rot) - P0*n_min_used/np.sqrt(Lambda)
        # Obtain index of minimal difference/distance
        index_max = np.where(abs(diff_max) == np.min(abs(diff_max)))[0]
        index_min = np.where(abs(diff_min) == np.min(abs(diff_min)))[0]
        # Calculate the rotationally shifted frequency (TAR approximation)
        ### in the inertial frame
        if rotation_frame == 'INERTIAL':
            f_min = (np.sqrt(Lambda[index_max]) / (P0*n_max_used) + m*omega_rot)[0]
            f_max = (np.sqrt(Lambda[index_min]) / (P0*n_min_used) + m*omega_rot)[0]
        ### in the co-rotating frame
        else:
            f_min = (np.sqrt(Lambda[index_max]) / (P0*n_max_used))[0]
            f_max = (np.sqrt(Lambda[index_min]) / (P0*n_min_used))[0]
    return f_min, f_max
################################################################################
################################################################################
# Function written by Jordan Van Beeck
################################################################################
def calculate_k(l,m,rossby):
  """
    Compute the mode classification parameter for gravity or Rossby modes from the corresponding azimuthal order (m) and spherical degree (l).
    Raises an error when l is smaller than m.
    ------- Parameters -------
    rossby: boolean
        parameter that needs to be set to True if Rossby mode k is calculated
    l, m: integer
        degree (l) and azimuthal order (m) of the modes
    ------- Returns -------
    k: integer
        mode classification parameter of the pulsation mode
  """
  if not rossby:
    # g-mode k
    if abs(l) >= abs(m):
      k = l - abs(m) # Lee & Saio (1997) (& GYRE source code --> see below)
      return k
    else:
      raise Exception(f'l is smaller than m, please revise your script/logic. The corresponding values were: (l,m) = ({l},{m})')
  else:
    # Rossby mode k
    if abs(l) >= abs(m):
      k = (-1)*(l - abs(m) + 1) # see GYRE source code: /gyre/src/build/gyre_r_tar_rot.f90 ; function r_tar_rot_t_ (Townsend & Teitler (2013))
      return k
    else:
      raise Exception(f'l is smaller than m, please revise your script/logic. The corresponding values were: (l,m) = ({l},{m})')

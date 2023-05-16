"""A few helpful functions to process MESA output."""
# from foam import functions_for_mesa as ffm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import multiprocessing, glob, h5py
from pathlib import Path
from functools import partial
from foam import support_functions as sf

################################################################################
def read_mesa_file(file_path, index_col=None):
    """
    Read in a mesa profile or history file and return 2 dictionaries.
    The first with the header info, and the second with the data.
    If the format is hdf5, assumes the attributes are the MESA header.
    If the format is not an hdf5 file, assumes the default MESA output format.
    This format is an ascii file delimited by whitespace with the following structure:
    line 1: header names
    line 2: header data
    line 3: blank
    line 4: main data names
    line >4:main data values

    ------- Parameters -------
    file_path: String
        The path to the MESA profile or history file to read in.
    ------- Returns -------
    header: dictionary
        A dictionary holding the header info of the MESA file.
    data: dictionary
        A dictionary holding the data columns of the MESA file as numpy arrays.
    """
    if h5py.is_hdf5(file_path):
        return sf.read_hdf5(file_path)

    else:   # assumes the default MESA output format
        header_df = pd.read_table(file_path, delim_whitespace=True, nrows=1, header=1)
        data_df = pd.read_table(file_path, delim_whitespace=True, skiprows=3, header=1, index_col=index_col)

        header={}
        for k in header_df.keys():
            header.update({k: header_df[k].to_numpy()[0]})
        data = {}
        for k in data_df.keys():
            data.update({k: data_df[k].to_numpy()})

        return header, data

################################################################################
def check_hydro_eq(profile_file, treshold_for_plot=5E-8):
    """
    Calculates the normalised differences between the terms on both sides of the hydrostatic equilibrium equation.
    Makes a plot if the values exceed a given treshold.
    This shows how well hydrostatic equilibrium has been fulfilled in the MESA model.
    ------- Parameters -------
    profile_file: String
        The path to the profile file to be checked.
    treshold_for_plot: float
        Make a plot if the max difference between the left and right hand side of the equation is greater than this number.
    """
    header, data = read_mesa_file(profile_file)
    # Compute the hydrostatic equilibrium quality factor q - See e.g. Aerts et al. (2010)
    lhs=np.delete(data['hyeq_lhs'], 0)  # remove the values at the surface, since these are 0
    rhs=np.delete(data['hyeq_rhs'], 0)

    norm = np.max(np.vstack([lhs, rhs]), axis=0)
    hyeq = np.abs(lhs - rhs) / np.abs(norm)

    # Make the plot if the treshold criterion is met
    if max(hyeq) > treshold_for_plot:
        # print the maximal deviation and profile name
        print(max(hyeq))
        print(profile_file)
        # make a semilog plot
        fig=plt.figure()
        ax = fig.add_subplot(111)
        ax.semilogy(np.delete(data['radius'], 0), hyeq, 'ko-') #also remove surface value, to have same amount of datapoints
        plt.show()
        plt.clf()
        plt.close('all')

################################################################################
def calculate_number_densities(hist_file):
    '''
    Calculate surface number densities for all isotopes in the MESA grid.
    All isotopes in the used nuclear network need to be be written in the history file,
    otherwise this will function give wrong numbers.
    ------- Parameters -------
    hist_file: String
        The path to the MESA history file.
    ------- Returns -------
    number_densities: dictionary
        Column keys specify the element (surf_X_per_N_tot), values are number densities of that element.
    '''
    header, data = read_mesa_file(hist_file)
    element_list = {}
    number_densities = {}
    inverse_average_atomic_mass = np.zeros(len(data[ list(data.keys())[0] ]))
    for column_name in data.keys():
        if '_per_Mass_tot' in column_name:
            element_list.update({column_name: data[column_name]})
            inverse_average_atomic_mass += data[column_name]

    average_atomic_mass = inverse_average_atomic_mass**(-1)
    for key in element_list.keys():
        number_densities.update({ key.replace('_per_Mass_tot', '_per_N_tot') : element_list[key]*average_atomic_mass})

    return number_densities

################################################################################
def convert_units(quantity, input, convertto='cgs'):
    '''
    Converts from solar units to cgs and vice versa.
    ------- Parameters -------
    input: list of float
        Numbers to convert.
    quantity: string
        The quantity you want to convert. Options are: {radius, mass, luminosity}.
    convertto: string
        The unit system to convert to. Options are: {cgs,solar} with default: 'cgs'.
    ------- Returns -------
    out: list of float
        converted numbers
    '''
    # conversion factors used in MESA see $MESA_DIR/const/public/const_def.f90
    to_cgs = {'mass': 1.9892E33,
              'luminosity': 3.8418E33,
              'radius': 6.9598E10,
              'cgrav': 6.67428E-8} # gravitational constant (g^-1 cm^3 s^-2)

    # convert to cgs
    if convertto == 'cgs':
        return input * to_cgs[quantity]
    # convert to solar units
    elif convertto == 'solar':
        return input / to_cgs[quantity]

################################################################################
def grid_extract_spectroscopy(mesa_profiles, output_file='gridSpectroscopy.hdf', parameters=['Z', 'M', 'logD', 'aov', 'fov', 'Xc'], nr_cpu=None):
    """
    Extract spectroscopic info and age for each globbed MESA profile and write them to 1 large file.
    ------- Parameters -------
    mesa_profiles: string
        String to glob to find all the relevant MESA profiles.
    output_file: string
        Name (or path) for the file containing all the pulsation frequencies of the grid.
    parameters: list of strings
        List of parameters varied in the computed grid, so these are taken from the
        name of the profile files, and included in the ouput file containing the info of the whole grid.
    nr_cpu: int
        Number of worker processes to use in multiprocessing. The default 'None' will use the number returned by os.cpu_count().
    """
    extract_func = partial(spectro_from_profiles, parameters=parameters)
    # Glob all the files, then iteratively send them to a pool of processors
    profiles = glob.iglob(mesa_profiles)
    with multiprocessing.Pool(nr_cpu) as p:
        spectro = p.imap(extract_func, profiles)

        # Generate the directory for the output file and write the file afterwards
        Path(Path(output_file).parent).mkdir(parents=True, exist_ok=True)
        # make a new list, so 'parameters' is not extended before passing it on to 'spectro_from_profiles'
        header_parameters = list(parameters)
        header_parameters.extend(['logTeff', 'logL', 'logg', 'age'])

        # Make list of lists, put it in a dataframe, and write to a file
        data = []
        for line in spectro:
            data.append(line)

    df = pd.DataFrame(data=data, columns=header_parameters)
    df.to_hdf(f'{output_file}', 'spectrogrid', format='table', mode='w')


################################################################################
def spectro_from_profiles(mesa_profile, parameters):
    """
    Extract spectroscopic info and age from a MESA profile and the model parameters from its filename.
    ------- Parameters -------
    mesa_profile: string
        path to the MESA profile
    parameters: list of strings
        List of input parameters varied in the computed grid, so these are read from the filename and included in returned line.

    ------- Returns -------
    line: string
        Line containing all the model- and spectroscopic parameters of the MESA profile.
    """
    param_dict = sf.get_param_from_filename(mesa_profile, parameters,values_as_float=True)
    prof_header, prof_data = read_mesa_file(mesa_profile)

    logL = np.log10(float(prof_header['photosphere_L']))
    logTeff = np.log10(float(prof_header['Teff']))
    logg = prof_data['log_g'][0]
    age=int(float(prof_header['star_age']))

    line=[]
    for p in parameters:
        line.append(param_dict[p])
    line.extend([logTeff, logL, logg, age])
    return line

################################################################################
def add_spectro_to_puls_grid(grid_frequencies, grid_spectroscopy, output_name='grid_spectro+freq.hdf', grid_parameters=['Z', 'M', 'logD', 'aov', 'fov', 'Xc']):
    """
    Combine the output files with the frequencies and spectroscopy of the grid in one new file,
    only keeping models that have entries in both the frequency and specto files.
    ------- Parameters -------
    grid_frequencies, grid_spectroscopy: string
        Paths to the files containing the model input parameters and corresponding frequency/spectroscopy of the model.
    output_name: string
        Name of the generated file containing the combined info.
    grid_parameters: list of string
        List of the model parameters to use for matching the entries in the freq/spectro file.
    """
    freq_df    = pd.read_hdf(grid_frequencies)
    spectro_df = pd.read_hdf(grid_spectroscopy)
    # Merge with spectro info first, freq info second. Only keeping rows that both dataFrames have in common based on the 'on' columns.
    df_merged  = pd.merge(spectro_df, freq_df, how='inner', on=grid_parameters)

    col = df_merged.pop("age") # Don't add the age in the combined file

    # take the column with rotation and place it as the first column, and its error as second column
    col = df_merged.pop("rot")
    df_merged.insert(0, col.name, col)
    col = df_merged.pop("rot_err")
    df_merged.insert(1, col.name, col)
    # write the merged dataFrame to a new file
    df_merged.to_hdf(f'{output_name}', 'puls_spectro_grid', format='table', mode='w')

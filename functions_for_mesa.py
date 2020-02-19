"""A few helpful functions to process MESA output."""
# from PyPulse import functions_for_mesa as ffm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

################################################################################
def read_mesa_file(file_path):
    """
    Read in a mesa profile or history file, and return 2 pandas dataframes.
    The first with the header info, and the second with the data.
    ------- Parameters -------
    file_path: String
        The path to the MESA profile or history file to read in.
    ------- Returns -------
    header: pandas dataframe
        A dataframe holding the header info of the MESA file.
    data: pandas dataframe
        A dataframe holding the data of the MESA file, the column used for indexing is
        the 'zone' for profiles, and the 'model_number' for the history file.
    """
    header = pd.read_table(file_path, delim_whitespace=True, nrows=1, header=1)
    data = pd.read_table(file_path, delim_whitespace=True, skiprows = 3, header=1, index_col=0)

    return header, data

################################################################################
def check_hydro_eq(profile_file, treshold_for_plot=5E-8):
    """
    Makes a plot showing the normalised differences between the terms on both sides of the hydrostatic equilibrium equation.
    This shows how well hydrostatic equilibrium has been fulfilled in your MESA model.
    ------- Parameters -------
    profile_file: String
        The path to the profile file to be used for the plot.
    treshold_for_plot: float
        Make a plot if the max difference between the left and right hand side of the equation is greater than this number.
    """
    header, data = read_mesa_file(profile_file)
    # Compute the hydrostatic equilibrium quality factor q - See e.g. Aerts et al. (2010)
    norm = np.max(np.vstack([data['hyeq_lhs'], data['hyeq_rhs']]), axis=0)
    hyeq = np.abs(data['hyeq_lhs'] - data['hyeq_rhs']) / np.abs(norm)
    hyeq_NaNs_removed = [h for h in hyeq if str(h) != 'nan']

    # Make the plot if the treshold criterion is met
    if max(hyeq_NaNs_removed) > treshold_for_plot:
        # print the maximal deviation and profile name
        print(max(hyeq_NaNs_removed))
        print(profile_file)
        # make a semilog plot
        fig=plt.figure()
        ax = fig.add_subplot(111)
        ax.semilogy(data['radius'], hyeq, 'ko-')
        plt.show()
    plt.close('all')

################################################################################
def plot_profile(profile_file, x_value, y_value, ax=None, label_size=16, colour='', linestyle='solid', alpha=1, legend=True, label=None):
    """
    Plot the requested quantities for the given MESA profile
    ------- Parameters -------
    profile_file: String
        The path to the profile file to be used for the plotting.
    x_value, y_value,: String
        The parameters of the profile plotted on the x and y axis.
        If x_value is mass or radius, it will be put in units relative to the total mass or radius
    ax: an axis object
        Axes object on which the plot will be made. If None: make figure and axis within this function.
    label_size, alpha: float
        The size of the labels in the figure, and transparency of the plot
    colour, linestyle, label: float
        Settings for the plot
    legend: boolean
        Flag to enable or disable a legend on the figure
    """
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    header, data = read_mesa_file(profile_file)
    # from "data", extract the columns
    y   = np.asarray(data[y_value])
    x   = np.asarray(data[x_value])
    if label == None:   #Set the label to be the name of the y variable
        label = y_value
    if x_value == 'radius' or x_value == 'mass':
        x = x/x[0] # normalized radius/mass coordinates
    # generate the plot, in which colour will not be specified
    if colour == '':
        ax.plot( x , y, label=label, linestyle=linestyle, alpha=alpha)
    # generate the plot, in which colour will be specified
    else:
        ax.plot( x , y, label=label, linestyle=linestyle, alpha=alpha, color=colour)
    if legend is True:
        ax.legend(loc='best', prop={'size': label_size})
    ax.set_xlabel(x_value, size=label_size)
    ax.set_ylabel(y_value, size=label_size)
################################################################################
################################################################################
def mesh_histogram(profile_file, x_value='radius', ax=None, label_size=16, colour='', linestyle='solid', alpha=1, legend=True, label=None, bins=200):
    """
    Make a histogram of the mesh points in the MESA profile
    ------- Parameters -------
    profile_file: String
        The path to the profile file to be used for the plotting.
    x_value: String
        The x value to use for the histogram
        If x_value is mass or radius, it will be put in units relative to the total mass or radius
    ax: an axis object
        Axes object on which the plot will be made. If None: make figure and axis within this function.
    label_size, alpha: float
        The size of the labels in the figure, and transparency of the plot
    colour, linestyle, label: float
        Settings for the plot
    legend: boolean
        Flag to enable or disable a legend on the figure.
    bins: int
        Number of bins used in the histogram
    """
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    header, data = read_mesa_file(profile_file)
    print(f'Total zones of {profile_file} : {header["num_zones"][0]}')

    # from "data", extract the columns
    x   = np.asarray(data[x_value])
    if label == None:   #Set the label to be the name of the y variable
        legend = False
    if x_value == 'radius' or x_value == 'mass':
        x = x/x[0] # normalized radius/mass coordinates
    # generate the plot, in which colour will not be specified
    if colour == '':
        ax.hist(x, bins=bins, histtype=u'step', label=label, alpha=alpha, linestyle=linestyle)
    # generate the plot, in which colour will be specified
    else:
        ax.hist(x, bins=bins, histtype=u'step', label=label, alpha=alpha, linestyle=linestyle, color=colour)
    # generate a legend if true
    if legend is True:
        ax.legend(loc='best', prop={'size': label_size})
    ax.set_xlabel(x_value, size=label_size)
    ax.set_ylabel('Meshpoints', size=label_size)

################################################################################
def plot_HRD(hist_file, ax=None, colour='blue', linestyle='solid', label='', label_size=16, Xc_marked=[0.1, 0.4, 0.7], Teff_logscale=True):
    """
    Makes an HRD plot from a provided history file
    ------- Parameters -------
    hist_file: String
        The path to the profile file to be used for the plot.
    ax: an axis object
        Axes object on which the plot will be made. If None: make figure and axis within this function.
    colour, linestyle, label: strings
        Specify the colour, linestyle and label of the plotted data
    label_size: float
        The size of the labels in the figure
    Xc_marked: list of floats
        Models with these Xc values are marked with red dots on the plot
    Teff_logscale: boolean
        Plot effective temperature in logscale (True), or not (False)
    """
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    header, data = read_mesa_file(hist_file)
    # from "data", extract the required columns as numpy arrays
    log_L     = np.asarray(data['log_L'])
    log_Teff  = np.asarray(data['log_Teff'])
    center_h1 = np.asarray(data['center_h1'])

    # Plot the x-axis in log scale
    if Teff_logscale:
        T = log_Teff
        ax.set_xlabel(r'log(T$_{eff}$)', size=label_size)
    # Plot the x-axis in linear scale
    else:
        T = 10**log_Teff
        ax.set_xlabel(r'T$_{eff}$ [K]', size=label_size)

    # Plot the HRD diagram (log_L vs. T)
    ax.plot( T, log_L, color = colour, linestyle = linestyle, label = label)
    ax.set_ylabel(r'log(L) [L$_{\odot}$]', size=label_size)
    ax.invert_xaxis()

    # Put specific marks on the HRD diagram
    k = 0
    for i in range(len(center_h1)-1, -1, -1):
        if center_h1[i] > Xc_marked[k]:
            ax.scatter( T[i] , log_L[i], marker='o', color = 'red', lw=2)
            k += 1
            if k >= len(Xc_marked):
                return
################################################################################
def calculate_number_densities(hist_file):
    '''
    Calculate surface number densities for all isotopes in the MESA grid.
    All isotopes in the used nuclear network need to be be written in the history file,
    or this will function give wrong numbers.
    ------- Parameters -------
    hist_file: String
        The path to the profile file to be used for the plot.
    ------- Returns -------
    number_densities: pandas dataframe
        Column keys specify the element (surf_X_per_N_tot), row keys the model number
        Values are number densities of that element
    '''
    header, data = read_mesa_file(hist_file)
    element_list = {}
    number_densities = {}
    inverse_average_atomic_mass = np.zeros(data.shape[0])
    for column_name in data.columns:
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

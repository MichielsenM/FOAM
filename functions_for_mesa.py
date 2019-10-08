# A few helpful functions to process MESA output.
# from PyPulse import functions_for_mesa as ffm
import matplotlib.pyplot as plt
import numpy as np
from . import read

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
    # Read the MESA profiles
    dic_prof  = read.read_multiple_mesa_files([profile_file], is_hist=False, is_prof=True)[0]
    # Extract data
    data      = dic_prof['prof']
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
def plot_profile(profile_file, x_value, y_value, ax=None, legend_size=18, colour='', linestyle='solid', alpha=1, legend=True):
    """
    Plot the requested quantities for the given MESA profile
    ------- Parameters -------
    profile_file: String
        The path to the profile file to be used for the plotting.
    x_value, y_value,: String
        The parameters of the profile plotted on the x and y axis.
        If x_value is mass or radius, it will be put in units relative to the total mass or radius
    ax: (optional) an axis object
        Axes object on which the plot will be made. If None: make figure and axis within this function.
    legend_size: float, optional
        The label size in the legend.
    """
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    dic_prof  = read.read_multiple_mesa_files([profile_file], is_hist=False, is_prof=True)[0] # create dictionary of MESA profile files
    # header    = dic_prof['header']
    
    # read the data
    data      = dic_prof['prof']

    # from "data", extract the columns
    y   = data[y_value]
    x   = data[x_value]
    if x_value == 'radius' or x_value == 'mass':
        x = x/x[0] # normalized radius/mass coordinates
    # generate the plot, in which colour will not be specified
    if colour == '':
        ax.plot( x , y, label=y_value, linestyle=linestyle, alpha=alpha)
    # generate the plot, in which colour will be specified
    else:
        ax.plot( x , y, label=y_value, linestyle=linestyle, alpha=alpha, color=colour)
    # generate a legend if true
    if legend is True:
        ax.legend(loc='best', prop={'size': legend_size})

################################################################################
def plot_HRD(hist_file, ax=None, colour='blue', ls='solid', label='', Xc_marked=[0.1, 0.4, 0.7], Teff_logscale=True):
    """
    Makes an HRD plot from a provided history file
    ------- Parameters -------
    hist_file: String
        The path to the profile file to be used for the plot.
    ax: (optional) an axis object
        Axes object on which the plot will be made. If None: make figure and axis within this function.
    colour, ls, label: (optional) strings
        Specify the colour, linestyle and label of the plotted data
    Xc_marked: (optional) list of floats
        Models with these Xc values are marked with red dots on the plot
    Teff_logscale: boolean
        Plot effective temperature in logscale (True), or not (False)
    """
    # if you don't receive a specified axes object, make one yourself
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    dic_hist  = read.read_multiple_mesa_files([hist_file], is_hist=True, is_prof=False)[0] # read the MESA profiles
    # extract the data
    data      = dic_hist['hist']
    # from "data", extract the required columns
    log_L     = data['log_L']
    log_Teff  = data['log_Teff']
    center_h1 = data['center_h1']

    # Plot the x-axis in log scale
    if Teff_logscale:
        T = log_Teff
        ax.set_xlabel(r'log(T$_{eff}$)') # grootheden onder een log kunnen geen eenheid hebben - JVB.
    # Plot the x-axis in linear scale
    else:
        T = 10**log_Teff
        ax.set_xlabel(r'T$_{eff}$ [K]')

    # Plot the HRD diagram (log_L vs. T)
    ax.plot( T, log_L, color = colour, linestyle = ls, label = label)
    ax.set_ylabel(r'log(L) [L$_{\odot}$]')
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
def convert_units(input, quantity, convertto='cgs'):
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
    # factors to convert
    to_cgs = {'mass': 1.9892E33,
              'luminosity': 3.8418E33,
              'radius': 6.9598E10}

    # convert to cgs
    if convertto == 'cgs':
        return input * to_cgs(quantity)
    # convert to solar units
    elif convertto == 'solar':
        return input / to_cgs(quantity)

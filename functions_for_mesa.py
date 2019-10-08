# A few helpful functions to process MESA output.
# from PyPulse import functions_for_mesa as ffm
import matplotlib.pyplot as plt
import numpy as np
from . import read

################################################################################
def check_hydro_eq(profile_file, treshold_for_plot=5E-8):
    """
    Makes a plot showing the normalised differences between both sides of the hydrostatic equilibrium equation.
    This shows how well hydrostatic equilibrium has been fulfilled.
    ------- Parameters -------
    profile_file: String
        The path to the profile file to be used for the plot.
    treshold_for_plot: float
        Make a plot if the max difference between the left and right hand side of the equation is above this.
    """
    dic_prof  = read.read_multiple_mesa_files([profile_file], is_hist=False, is_prof=True)[0]
    data      = dic_prof['prof']
    # Compute the hydrostatic equilibrium quality factor q
    norm = np.max(np.vstack([data['hyeq_lhs'], data['hyeq_rhs']]), axis=0)
    hyeq = np.abs(data['hyeq_lhs'] - data['hyeq_rhs']) / np.abs(norm)
    hyeq_NaNs_removed = [h for h in hyeq if str(h) != 'nan']

    if max(hyeq_NaNs_removed) > treshold_for_plot:
        print(max(hyeq_NaNs_removed))
        print(profile_file)
        # make a plot if difference becomes larger than specified value
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
        The path to the profile file to be used for the plot.
    x_value, y_value,: String
        The parameters from the profile to put on the x and y axis.
        If x_value is mass or radius, it will be put in units relative to the total mass or radius
    ax: (optional) an axis object
        Axes to make the plot on, make figure and axis in this function if None
    legend_size: float, optional
        The size of the labels in the legend.
    """
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    dic_prof  = read.read_multiple_mesa_files([profile_file], is_hist=False, is_prof=True)[0]
    # header    = dic_prof['header']
    data      = dic_prof['prof']

    # from "data", extract the columns
    y   = data[y_value]
    x   = data[x_value]
    if x_value == 'radius' or x_value == 'mass':
        x = x/x[0]

    if colour == '':
        ax.plot( x , y, label=y_value, linestyle=linestyle, alpha=alpha)
    else:
        ax.plot( x , y, label=y_value, linestyle=linestyle, alpha=alpha, color=colour)

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
        Axes to make the plot on, make figure and axis in this function if None
    colour, ls, label: (optional) strings
        Specify the colour, linestyle and label of the plotted data
    Xc_marked: (optional) list of floats
        Models with these Xc values are marked with red dots on the plot
    Teff_logscale: boolean
        Plot effective temperature in logscale, or not
    """
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    dic_hist  = read.read_multiple_mesa_files([hist_file], is_hist=True, is_prof=False)[0]
    data      = dic_hist['hist']
    # from "data", extract the columns
    log_L     = data['log_L']
    log_Teff  = data['log_Teff']
    center_h1 = data['center_h1']

    if Teff_logscale:
        T = log_Teff
        ax.set_xlabel(r'log(T$_{eff}$) [K]')
    else:
        T = 10**log_Teff
        ax.set_xlabel(r'T$_{eff}$ [K]')

    ax.plot( T, log_L, color = colour, linestyle = ls, label = label)
    ax.set_ylabel(r'log(L) [L$_{\odot}$]')
    ax.invert_xaxis()

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
    Converts from solar to cgs and vice versa.
    ------- Parameters -------
    input: list of float
        numbers to convert
    quantity: string
        which quantity you want to convert. Options are radius, mass, luminosity
    convertto: string
        the unit system to convert to, default is cgs but can choose solar.
    ------- Returns -------
    out: list of float
        converted numbers
    '''

    to_cgs = {'mass': 1.9892E33,
              'luminosity': 3.8418E33,
              'radius': 6.9598E10}

    if convertto == 'cgs':
        return input * to_cgs(quantity)

    elif convertto == 'solar':
        return input / to_cgs(quantity)

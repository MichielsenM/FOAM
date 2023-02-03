"""Helpful functions in general. Making figures, reading HDF5, processing strings."""
# from foam import support_functions as sf
import h5py, re
import pandas as pd
from pathlib import Path
import logging
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

logger = logging.getLogger('logger.sf')
################################################################################
def split_line(line, sep) :
    """
    Splits a string in 2 parts.

    ------- Parameters -------
    line: string
        String to split in 2.
    sep: string
        Separator where the string has to be split around.

    ------- Returns -------
    head: string
        Part 1 of the string before the separator.
    tail: string
        Part 2 of the string after the separator.
    """
    head, sep_, tail = line.partition(sep)
    assert sep_ == sep
    return head, tail

################################################################################
def substring(line, sep_first, sep_second) :
    """
    Get part of a string between 2 specified separators.
    If second separator is not found, return everyting after first separator.

    ------- Parameters -------
    line: string
        String to get substring from.
    sep_first: string
        First separator after which the returned substring should start.
    sep_first: string
        Second separator at which the returned substring should end.

    ------- Returns -------
    head: string
        Part of the string between the 2 separators.
    """
    head, tail = split_line(line, sep = sep_first)
    if sep_second not in tail:
        return tail
    head, tail = split_line(tail, sep =sep_second)
    return head
################################################################################
def get_param_from_filename(file_path, parameters):
    """
    Get parameters from filename

    ------- Parameters -------
    file_path : string
        Full path to the file
    parameters: list of Strings
        Names of parameters to extract from filename

    ------- Returns -------
    param_dict: Dictionary
        Keys are strings describing the parameter, values are strings giving corresponding parameter values
    """

    param_dict = {}
    for parameter in parameters:
        try:
            p = substring(Path(file_path).stem, parameter, '_')
            param_dict[parameter] = p
        except:
            param_dict[parameter] = '0'
            logger.info(f'In get_param_from_filename: parameter "{parameter}" not found in \'{file_path}\', value set to zero')

    return param_dict

################################################################################
def read_hdf5(filename):
    """
    Read a HDF5-format file (e.g. GYRE)

    ------- Parameters -------
    filename : string
        Input file

    ------- Returns -------
    attributes: dictionary
        Dictionary containing the attributes of the file.
    data: dictionary
        Dictionary containing the data from the file as numpy arrays.
    """
    # Open the file
    with h5py.File(filename, 'r') as file:
        # Read attributes
        attributes = dict(zip(file.attrs.keys(),file.attrs.values()))
        # Read datasets
        data = {}
        for k in file.keys() :
            data[k] = file[k][...]
    return attributes, data
################################################################################
def make_multipanel_plot(nr_panels=1, xlabel='', ylabels=[''], keys=None, title='', label_size=22, xlim=[],
                        left_space=0.1, bottom_space=0.085, right_space=0.978, top_space=0.97, h_space=0.12, figure_size = [12,8]):
    """
    Make a plot

    ------- Parameters -------
    nr_panels: int
        The number of panels to add to the plot
    xlabel, ylabels: string and list of strings
        Names for x and y axes.
    keys:
        The keys corresponding to the dictionary entries of the axes
    title: string
        Title of the plot, no title if argument not given.
    left_space, bottom_space, right_space, top_space: float, optional
        The space that needs to be left open around the figure.
    h_space: float
        The size of the space in between both panels.
    figure_size: list of 2 floats
        Specify the dimensions of the figure in inch.
    label_size: float
        The size of the labels on the axes, and slightly smaller on the ticks
    xlim: list
        Lists of length 2, specifying the lower and upper limits of the x axe.
        Matplotlib sets the limits automatically if the arguments are not given.

    ------- Returns -------
    ax_dict: Dictionary of the axes
    fig: The figure
    """
    if keys==None:   # make keys integers
        keys = range(nr_panels)

    fig=plt.figure(figsize=(figure_size[0], figure_size[1]))
    gs=GridSpec(nr_panels,1) # multiple rows, 1 column
    ax_dict = {}

    for i in range(0, nr_panels):
        if i==0:
            ax = fig.add_subplot(gs[i:i+1,0])
            if len(xlim) ==2:
                ax.set_xlim(xlim[0],xlim[1])

        else:
            ax = fig.add_subplot(gs[i:i+1,0], sharex=ax_dict[keys[0]])

        ax_dict.update({keys[i]:ax})

        ax.set_ylabel(ylabels[i], size=label_size)
        ax.tick_params(labelsize=label_size-2)

        if i == nr_panels-1:
            ax.set_xlabel(xlabel, size=label_size)
        else:
            plt.setp(ax.get_xticklabels(), visible=False)

    ax_dict[keys[0]].set_title(title)
    plt.subplots_adjust(hspace=h_space, left=left_space, right=right_space, top=top_space, bottom=bottom_space)

    return ax_dict, fig

################################################################################
def sign(x):
    """
    Returns the sign of a number as a string

    ------- Parameters -------
    x: float or int

    ------- Returns -------
    A string representing the sign of the number
    """
    if abs(x) == x:
        return '+'
    else:
        return '-'
################################################################################
def get_subgrid_dataframe(file_to_read, fixed_params=None):
    """
    Read a tsv file containing the grid information as a pandas dataframe.
    Parameters can be fixed to certain values to fiter out entries with other values of that parameter.
    ------- Parameters -------
    file_to_read: string
        path to the file to read
    fixed_params: dictionary
        keys are parameters to fix to the value specified in the dictionary

    ------- Returns -------
    df: pandas dataframe
    """
    df = pd.read_table(file_to_read, delim_whitespace=True, header=0)

    if fixed_params is not None:
        for param in fixed_params.keys():
            indices_to_drop = df[df[param] != fixed_params[param] ].index
            df.drop(indices_to_drop, inplace = True)
        df.reset_index(drop=True, inplace=True)

    return df

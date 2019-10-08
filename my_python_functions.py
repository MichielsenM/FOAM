# Helpful functions in general. Making figures, reading HDF5, processing strings.
# from PyPulse import my_python_functions as mypy
import h5py, re
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from . import read

import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('logger')
# logger.addHandler(logging.StreamHandler())
logger.setLevel(logging.DEBUG)

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
        Strings specifying the given parameters
    """

    param_dict = {}
    path, filename = file_path.rsplit('/',1)
    filename = filename[:filename.rfind('.')] # remove file extension

    for parameter in parameters:
        try:
            head, tail = split_line(filename, parameter)
            # p, tail =  re.split('\\.|_', tail, 1) # alternative way
            p = substring(filename, parameter, '_')
            param_dict[parameter] = p
        except:
            param_dict[parameter] = '0'
            logger.info(f'In get_param_from_filename: parameter "{parameter}" not found in filename, value set to zero')

    return param_dict

################################################################################
def read_hdf5(filename):
    """
    Read a HDF5-format file (e.g. GYRE)

    ------- Parameters -------
    filename : string
        Input file

    ------- Returns -------
    data: dictionary
        Dictionary containing the information from the file.
    """
    # Open the file
    with h5py.File(filename, 'r') as file:
        # Read attributes
        data = dict(zip(file.attrs.keys(),file.attrs.values()))
        # Read datasets
        for k in file.keys() :
            data[k] = file[k][...]

    return data
################################################################################
def make_multipanel_plot(nr_panels=1, xlabel='', ylabels=[''], keys=None, title='', left_space=0.1, bottom_space=0.085, right_space=0.978, top_space=0.97, h_space=0.12, figure_size = [12,8], label_size=22, xlim=[]):
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

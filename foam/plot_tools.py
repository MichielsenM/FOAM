""" Plotting functionality """
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator
from pathlib import Path
from foam import functions_for_mesa as ffm
################################################################################
def make_multipanel_plot(nr_panels=1, xlabel='', ylabels=[''], keys=None, title='', label_size=22, xlim=[],
                        left_space=0.1, bottom_space=0.085, right_space=0.978, top_space=0.97, h_space=0.12, figure_size = [12,8]):
    """
    Make a multipanel figure for plots.

    Parameters
    ----------
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

    Returns
    ----------
    ax_dict: dict 
        Dictionary of the axes
    fig: Figure
        The figure
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
def corner_plot(merit_values_file, merit_values_file_error_ellipse, fig_title, observations_file, label_size=20, fig_output_dir='figures_correlation/',
                      percentile_to_show=0.5, logg_or_logL='logL', mark_best_model= False, n_sigma_box=3, grid_parameters = None,
                      axis_labels_dict = {'rot': r'$\Omega_{\mathrm{rot}}$ [d$^{-1}$]' ,'M': r'M$_{\rm ini}$', 'Z': r'Z$_{\rm ini}$',
                      'logD':r'log(D$_{\rm env}$)', 'aov':r'$\alpha_{\rm CBM}$','fov':r'f$_{\rm CBM}$','Xc':r'$\rm X_c$'} ):
    """
    Make a plot of all variables vs each other variable, showing the MLE values as colorscale.
    A kiel/HR diagram is made, depending on if logg_obs or logL_obs is passed as a parameter.
    The subplots on the diagonal show the distribution of that variable.
    The list of variables is retrieved from columns of the merit_values_file,
    where the first column is 'meritValue', which are the MLE values.
    The resulting figure is saved afterwards in the specified location.

    Parameters
    ----------
    merit_values_file, merit_values_file_error_ellipse: string
        Path to the hdf5 files with the merit function values and parameters of the models in the grid,
        and of just the models in the error ellipse.
    observations_file: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    fig_title: string
        Title of the figure and name of the saved png.
    label_size: int
        Size of the axis labels.
    fig_output_dir: string
        Output directory for the figures.
    percentile_to_show: float
        Percentile of models to show in the plots.
    logg_or_logL: string
        String 'logg' or 'logL' indicating whether log of surface gravity (g) or luminosity (L) is plot.
    mark_best_model: boolean
        Indicate the best model with a marker
    grid_parameters: list of string
        List of the parameters in the theoretical grid.
    axis_labels_dict: dictionary
        Keys are grid parameters, values are strings how those values should be shown on the axis labels
    """
    # Define custom colormap
    cdict = {'red':((0.0, 1.0, 1.0),
                    (0.5, 1.0, 1.0),
                    (0.75, 1.0, 1.0),
                    (1.0, 0.75, 0.75)),
             'green': ((0.0, 1.0, 1.0),
                       (0.25, 0.5, 0.5),
                       (0.5, 0.0, 0.0),
                       (1.0, 0.0, 0.0)),
             'blue':  ((0.0, 0.0, 0.0),
                       (0.5, 0.0, 0.0),
                       (1.0, 1.0, 1.0))
            }
    CustomCMap = LinearSegmentedColormap('CustomMap', cdict)
    # theoretical models within the error ellipse
    dataframe_theory_error_ellipse = pd.read_hdf(merit_values_file_error_ellipse, delim_whitespace=True, header=0)
    dataframe_theory_error_ellipse = dataframe_theory_error_ellipse.sort_values('meritValue', ascending=False)    # Order from high to low, to plot lowest values last

    # theoretical models
    dataframe_theory = pd.read_hdf(merit_values_file)
    dataframe_theory = dataframe_theory.sort_values('meritValue', ascending=False)    # Order from high to low, to plot lowest values last
    dataframe_theory = dataframe_theory.iloc[int(dataframe_theory.shape[0]*(1-percentile_to_show)):] # only plot the given percentage lowest meritValues

    if dataframe_theory.iloc[0]['rot'] == dataframe_theory.iloc[1]['rot'] == dataframe_theory.iloc[2]['rot'] == dataframe_theory.iloc[-1]['rot']: # rotation is fixed, don't plot it
        df_error_ellipse = dataframe_theory_error_ellipse.filter(['meritValue']+grid_parameters ) # make new dataframe with only needed info
        df = dataframe_theory.filter(['meritValue']+grid_parameters )       # make new dataframe with only needed info
        # Remove models in the error ellipse from the regular dataframe.
        df = pd.merge(df,df_error_ellipse, indicator=True, how='outer', on=grid_parameters, suffixes=[None, '_remove']).query('_merge=="left_only"').drop(['meritValue_remove', '_merge'], axis=1)

    else: # rotation was varied, include it in the plots
        df_error_ellipse = dataframe_theory_error_ellipse.filter(['meritValue']+['rot']+grid_parameters ) # make new dataframe with only needed info
        df = dataframe_theory.filter(['meritValue']+['rot']+grid_parameters )       # make new dataframe with only needed info
        # Remove models in the error ellipse from the regular dataframe.
        df = pd.merge(df,df_error_ellipse, indicator=True, how='outer', on=grid_parameters, suffixes=[None, '_remove']).query('_merge=="left_only"').drop(['meritValue_remove', 'rot_remove', '_merge'], axis=1)

    ax_dict={}  # dictionary of dictionaries, holding the subplots of the figure, keys indicate position (row, column) of the subplot
    nr_params = len(df.columns)-1
    for i in range(nr_params):
        ax_dict.update({i:{}})

    fig=plt.figure(figsize=(10,8))
    gs=GridSpec(nr_params,nr_params) # multiple rows and columns

    if mark_best_model: min_index = df_error_ellipse['meritValue'].idxmin(axis='index', skipna=True)    # get the best model according to the point estimator

    for ix in range(0, nr_params):
        for iy in range(0, nr_params-ix):
            if iy==0:
                share_x = None
            else:
                share_x = ax_dict[0][ix]
            if (ix==0) or (iy+ix == nr_params-1):
                share_y = None
            else:
                share_y = ax_dict[iy][0]

            # create subplots and add them to the dictionary
            ax = fig.add_subplot(gs[nr_params-iy-1:nr_params-iy,ix:ix+1], sharex=share_x, sharey=share_y)
            ax_dict[iy].update({ix:ax})

            # manage visibility and size of the labels and ticks
            ax.tick_params(labelsize=label_size-4)
            if ix == 0:
                ax.set_ylabel(axis_labels_dict[df.columns[iy+1]], size=label_size)
                if (iy == nr_params-1):
                    plt.setp(ax.get_yticklabels(), visible=False)
            else:
                plt.setp(ax.get_yticklabels(), visible=False)
            if iy == 0:
                ax.set_xlabel(axis_labels_dict[df.columns[nr_params-ix]], size=label_size)
                ax.tick_params(axis='x', rotation=45)
            else:
                plt.setp(ax.get_xticklabels(), visible=False)

            if (iy+ix == nr_params-1):  # make distribution plots on the diagonal subplots
                values = sorted(np.unique(df.iloc[:,nr_params-ix]))
                # determine edges of the bins for the histogram distribution plots
                if df.columns[nr_params-ix] == 'rot':
                    domain = (values[0], values[-1])
                    ax.hist( df_error_ellipse.iloc[:,nr_params-ix], bins=25, range=domain, density=False, cumulative=False, histtype='step' )

                else:
                    if len(values) > 1:
                        bin_half_width = (values[0]+values[1])/2-values[0]
                    else:
                        bin_half_width = 1E-3
                    bin_edges = [values[0]-bin_half_width]
                    for i in range(len(values)-1):
                        bin_edges.extend([(values[i]+values[i+1])/2])
                    bin_edges.extend([values[-1]+bin_half_width])
                    ax.hist( df_error_ellipse.iloc[:,nr_params-ix], bins=bin_edges, density=False, cumulative=False, histtype='step' )

                ax.tick_params(axis='y',left=False)
                continue

            im = ax.scatter(df.iloc[:,nr_params-ix], df.iloc[:,iy+1], c=np.log10(df.iloc[:,0]), cmap='Greys_r')
            im = ax.scatter(df_error_ellipse.iloc[:,nr_params-ix], df_error_ellipse.iloc[:,iy+1], c=np.log10(dataframe_theory_error_ellipse['meritValue']), cmap=CustomCMap)
            if mark_best_model: ax.scatter(df_error_ellipse.loc[min_index][nr_params-ix], df.loc[min_index][iy+1], color='white', marker = 'x')
            # Adjust x an y limits of subplots
            limit_adjust = (max(df.iloc[:,iy+1]) - min(df.iloc[:,iy+1]))*0.08
            if limit_adjust == 0: limit_adjust=0.1
            ax.set_ylim( min(df.iloc[:,iy+1])-limit_adjust,  max(df.iloc[:,iy+1])+limit_adjust  )
            limit_adjust = (max(df.iloc[:,nr_params-ix]) - min(df.iloc[:,nr_params-ix])) *0.08
            if limit_adjust == 0: limit_adjust=0.1
            ax.set_xlim( min(df.iloc[:,nr_params-ix])-limit_adjust, max(df.iloc[:,nr_params-ix])+limit_adjust )

    fig.align_labels()
    # add subplot in top right for Kiel or HRD
    ax_hrd = fig.add_axes([0.508, 0.65, 0.33, 0.33]) # X, Y, width, height

    ax_hrd.set_xlabel(r'log(T$_{\mathrm{eff}}$ [K])', size=label_size)
    ax_hrd.tick_params(labelsize=label_size-4)
    ax_hrd.invert_xaxis()

    # Observations
    if n_sigma_box != None:
        obs_dataframe  = pd.read_table(observations_file, delim_whitespace=True, header=0, index_col='index')
        if (('logL' in obs_dataframe.columns) or ('logg' in obs_dataframe.columns)) and ('Teff' in obs_dataframe.columns)  :
            if 'logL' not in obs_dataframe.columns:
                logg_or_logL = 'logg'

            # Observed spectroscopic error bar, only added if observational constraints were provided.
            # To add the 1 and n-sigma spectro error boxes, calculate their width (so 2 and 2*n sigma wide)
            width_logTeff_sigma= np.log10(obs_dataframe['Teff'][0]+obs_dataframe['Teff_err'][0]) - np.log10(obs_dataframe['Teff'][0]-obs_dataframe['Teff_err'][0])
            width_logTeff_nsigma= np.log10(obs_dataframe['Teff'][0]+n_sigma_box*obs_dataframe['Teff_err'][0]) - np.log10(obs_dataframe['Teff'][0]-n_sigma_box*obs_dataframe['Teff_err'][0])
            errorbox_1s = patches.Rectangle((np.log10(obs_dataframe['Teff'][0]-obs_dataframe['Teff_err'][0]),obs_dataframe[logg_or_logL][0]-obs_dataframe[f'{logg_or_logL}_err'][0]),
                        width_logTeff_sigma, 2*obs_dataframe[f'{logg_or_logL}_err'][0],linewidth=1.7,edgecolor='cyan',facecolor='none', zorder=2.1)
            errorbox_ns = patches.Rectangle((np.log10(obs_dataframe['Teff'][0]-n_sigma_box*obs_dataframe['Teff_err'][0]), obs_dataframe[logg_or_logL][0]-n_sigma_box*obs_dataframe[f'{logg_or_logL}_err'][0]),
                        width_logTeff_nsigma, 2*n_sigma_box*obs_dataframe[f'{logg_or_logL}_err'][0],linewidth=1.7,edgecolor='cyan',facecolor='none', zorder=2.1)
            ax_hrd.add_patch(errorbox_1s)
            ax_hrd.add_patch(errorbox_ns)

    if logg_or_logL=='logg': ax_hrd.invert_yaxis()

    im = ax_hrd.scatter(dataframe_theory['logTeff'], dataframe_theory[logg_or_logL], c=np.log10(dataframe_theory['meritValue']), cmap='Greys_r')
    im_error_ellipse = ax_hrd.scatter(dataframe_theory_error_ellipse['logTeff'], dataframe_theory_error_ellipse[logg_or_logL], c=np.log10(dataframe_theory_error_ellipse['meritValue']), cmap=CustomCMap)
    ax_hrd.set_ylabel(f'{logg_or_logL[:-1]} {logg_or_logL[-1]}')
    if logg_or_logL == 'logL':
        ax_hrd.set_ylabel(r'log(L/L$_{\odot}$)', size=label_size)
    elif logg_or_logL == 'logg':
        ax_hrd.set_ylabel(r'$\log\,g$ [dex]', size=label_size)

    ax_hrd.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax_hrd.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax_hrd.tick_params(which='major', length=6)
    ax_hrd.tick_params(which='minor', length=4)

    if mark_best_model: ax_hrd.scatter(dataframe_theory_error_ellipse['logTeff'][min_index], dataframe_theory_error_ellipse[logg_or_logL][min_index], marker='x', color='white')

    # Add color bar
    cax = fig.add_axes([0.856, 0.565, 0.04, 0.415]) # X, Y, width, height
    cbar= fig.colorbar(im, cax=cax, orientation='vertical')
    cax2 = fig.add_axes([0.856, 0.137, 0.04, 0.415]) # X, Y, width, height
    cbar2= fig.colorbar(im_error_ellipse, cax=cax2, orientation='vertical', )

    if dataframe_theory_error_ellipse.shape[0]==1: # To prevent messing up colours due to automatic rescaling of colorbar
        im_error_ellipse.set_clim(np.log10(dataframe_theory_error_ellipse['meritValue']), np.log10(dataframe_theory_error_ellipse['meritValue'])*1.1)

    if '_MD_' in fig_title:
        cbar.set_label('log(MD)', rotation=90, size=label_size)
        cbar2.set_label('log(MD)', rotation=90, size=label_size)
    elif '_CS_' in fig_title:
        cbar.set_label(r'log($\chi^2$)', rotation=90, size=label_size)
        cbar2.set_label(r'log($\chi^2$)', rotation=90, size=label_size)
    else:
        cbar.set_label('log(merit function value)', rotation=90)
    cbar.ax.tick_params(labelsize=label_size-4)
    cbar2.ax.tick_params(labelsize=label_size-4)
    fig.subplots_adjust(left=0.114, right=0.835, bottom=0.137, top=0.99)

    # fig.suptitle(fig_title, horizontalalignment='left', size=20, x=0.28)
    Path(fig_output_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig(f'{fig_output_dir}{fig_title}.png', dpi=400)
    plt.clf()
    plt.close(fig)
################################################################################
def plot_mesa_file(profile_file, x_value, y_value, ax=None, label_size=16, colour='', linestyle='solid', alpha=1, legend=True, label=None):
    """
    Plot the requested quantities for the given MESA profile or history file.

    Parameters
    ----------
    profile_file: string
        The path to the profile file to be used for the plotting.
    x_value, y_value,: string
        The parameters of the profile plotted on the x and y axis.
        If x_value is mass or radius, it will be put in units relative to the total mass or radius
    ax: Axes
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

    header, data = ffm.read_mesa_file(profile_file)
    # from "data", extract the columns
    y   = np.asarray(data[y_value])
    x   = np.asarray(data[x_value])
    if label == None:   #Set the label to be the name of the y variable
        label = y_value
    if x_value == 'radius' or x_value == 'mass':
        x = x/x[0] # normalized radius/mass coordinates
        ax.set_xlim(0,1)
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
def plot_mesh_histogram(profile_file, x_value='radius', ax=None, label_size=16, colour='', linestyle='solid', alpha=1, legend=True, label=None, bins=200):
    """
    Make a histogram of the mesh points in the MESA profile.

    Parameters
    ----------
    profile_file: string
        The path to the profile file to be used for the plotting.
    x_value: string
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

    header, data = ffm.read_mesa_file(profile_file)
    print(f'Total zones of {profile_file} : {header["num_zones"]}')

    # from "data", extract the columns
    x   = np.asarray(data[x_value])
    if label == None:   #Set the label to be the name of the y variable
        legend = False
    if x_value == 'radius' or x_value == 'mass':
        x = x/x[0] # normalized radius/mass coordinates
        ax.set_xlim(0,1)
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
def plot_hrd(hist_file, ax=None, colour='blue', linestyle='solid', label='', label_size=16,
    Xc_marked=None, Teff_logscale=True, start_track_from_Xc=None, diagram='HRD'):
    """
    Makes an HRD plot from a provided MESA history file.

    Parameters
    ----------
    hist_file: string
        The path to the profile file to be used for the plot.
    ax: Axes
        Axes object on which the plot will be made. If None: make figure and axis within this function.
    colour, linestyle, label: strings
        Specify the colour, linestyle and label of the plotted data.
    label_size: float
        The size of the labels in the figure.
    Xc_marked: list of floats
        Models with these Xc values are marked with red dots on the plot (listed in increasing value).
    Teff_logscale: boolean
        Plot effective temperature in logscale (True), or not (False).
    start_track_from_Xc: float
        Only start plotting the track if Xc drops below this value (e.g. to not plot the initial relaxation loop).
    diagram: string
        Type of diagram that is plotted. Options are HRD (logL vs logTeff), sHRD (log(Teff^4/g) vs logTeff) or kiel (logg vs logTeff).
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)

    header, data = ffm.read_mesa_file(hist_file)

    # From "data", extract the required columns as numpy arrays
    log_L     = np.asarray(data['log_L'])
    log_Teff  = np.asarray(data['log_Teff'])
    log_g     = np.asarray(data['log_g'])
    center_h1 = np.asarray(data['center_h1'])
    # Plot the x-axis in log scale
    if Teff_logscale:
        T = log_Teff
        ax.set_xlabel(r'log(T$_{\mathrm{eff}}$)', size=label_size)
    # Plot the x-axis in linear scale
    else:
        T = 10**log_Teff
        ax.set_xlabel(r'T$_{\mathrm{eff}}$ [K]', size=label_size)

    # Plot HRD
    if diagram == 'HRD':
        y_axis = log_L
        ax.set_ylabel(r'log(L/L$_{\odot}$)', size=label_size)
    # Plot sHRD (log_Teff^4/log_g vs log_Teff)
    elif diagram == 'sHRD':
        log_Lsun = 10.61
        y_axis = 4*log_Teff - log_g - log_Lsun
        ax.set_ylabel(r'$\log \left(\frac{{T_{\mathrm{eff}}}^4}{g}\right) \ (\mathscr{L}_\odot)$', size=label_size)
    # Plot Kiel diagram (log_g vs log_Teff)
    elif diagram == 'kiel':
        y_axis = log_g
        ax.set_ylabel(r'log g [dex]', size=label_size)

    # Start plotting from Xc value
    if start_track_from_Xc!= None:
        for i in range(len(center_h1)):
            if center_h1[i] < start_track_from_Xc:
                T = T[i:]
                y_axis = y_axis[i:]
                break

    # Plot the HRD diagram (log_L vs. T)
    ax.plot(T, y_axis, color = colour, linestyle = linestyle, label = label)

    # Put specific marks on the HRD diagram
    if Xc_marked is None:
        return
    k = 0
    for i in range(len(center_h1)-1, -1, -1):
        if center_h1[i] > Xc_marked[k]:
            ax.scatter( T[i] , log_L[i], marker='o', color = 'red', lw=2)
            k += 1
            if k >= len(Xc_marked):
                return

################################################################################
def plot_khd(hist_file, ax=None, number_mix_zones=8, xaxis='model_number'):
    """
    Makes a Kippenhahn plot from a provided MESA history file.

    Parameters
    ----------
    hist_file: string
        The path to the history file to be used for the plot.
    ax: Axes
        Axes object on which the plot will be made. If None: make figure and axis within this function.
    number_mix_zones: int
        Number of mixing zones included in the mesa history file.
    xaxis: string
        Quantity to put on the x-axis of the plot (e.g. model_number or star_age).
    """
    if ax is None:
        fig = plt.figure( figsize=(10, 4))
        ax = fig.add_subplot(111)

    _, data = ffm.read_mesa_file(hist_file)

    x_values=data[xaxis]
    m_star = data['star_mass']
    m_ini  = m_star[0]

    for j in range(number_mix_zones):
        colours = {'-1':'w', '0':'w', '1':'lightgrey', '2':'b', '7':'g', '3':'cornflowerblue', '8':'red'}
        if j == number_mix_zones-1:
            ax.vlines(x_values, 0, data[f'mix_qtop_{number_mix_zones-j}']*m_star/m_ini, color=[colours[str(x)] for x in data[f'mix_type_{number_mix_zones-j}']])
        else:
            ax.vlines(x_values, data[f'mix_qtop_{number_mix_zones-1-j}']*m_star/m_ini, data[f'mix_qtop_{number_mix_zones-j}']*m_star/m_ini, color=[colours[str(x)] for x in data[f'mix_type_{number_mix_zones-j}']])

    ax.plot(x_values, m_star / m_ini, lw=1, color='black', label=f'{m_ini:.1f} $M_\odot$')

    ax.set_xlim(min(x_values)*0.99, max(x_values)*1.01)
    ax.set_ylim(0, 1.02)

    ax.plot([], [], lw=10, color='lightgrey', label=r'Convective') # Only to make it appear in the Legend
    ax.plot([], [], lw=10, color='cornflowerblue', label=r'CBM')   # Only to make it appear in the Legend
    ax.legend(bbox_to_anchor=(0.15, 0.97), fontsize=10, frameon=False, fancybox=False, shadow=False, borderpad=False)
    ax.set_ylabel(r'Relative Mass $\, m/M_\star$')
    ax.set_xlabel(xaxis)
    plt.tight_layout()

    return
################################################################################

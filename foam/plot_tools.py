import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import AutoMinorLocator
from pathlib import Path

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
def corner_plot(merit_values_file, merit_values_file_error_ellips, fig_title, observations_file, label_size=20, fig_outputDir='figures_correlation/',
                      percentile_to_show=0.5, logg_or_logL='logL', mark_best_model= False, n_sigma_spectrobox=3):
    """
    Make a plot of all variables vs each other variable, showing the MLE values as colorscale.
    A kiel/HR diagram is made, depending on if logg_obs or logL_obs is passed as a parameter.
    The subplots on the diagonal show the distribution of that variable.
    The list of variables is retrieved from columns of the merit_values_file,
    where the first column is 'meritValue', which are the MLE values.
    The resulting figure is saved afterwards in the specified location.
    ------- Parameters -------
    merit_values_file, merit_values_file_error_ellips: string
        Path to the tsv files with the merit function values and parameters of the models in the grid,
        and of just the models in the error ellips.
    observations_file: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    fig_title: string
        Title of the figure and name of the saved png.
    label_size: int
        Size of the axis labels.
    fig_outputDir: string
        Output directory for the figures.
    percentile_to_show: float
        Percentile of models to show in the plots.
    logg_or_logL: string
        String 'logg' or 'logL' indicating wheter log of surface gravity (g) or luminosity (L) is plot.
    mark_best_model: boolean
        Indicate the best model with a marker
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
    # theoretical models within the error ellips
    df_Theo_EE = pd.read_table(merit_values_file_error_ellips, delim_whitespace=True, header=0)
    df_Theo_EE = df_Theo_EE.sort_values('meritValue', ascending=False)    # Order from high to low, to plot lowest values last

    # theoretical models
    df_Theo = pd.read_hdf(merit_values_file)
    df_Theo = df_Theo.sort_values('meritValue', ascending=False)    # Order from high to low, to plot lowest values last
    df_Theo = df_Theo.iloc[int(df_Theo.shape[0]*(1-percentile_to_show)):] # only plot the given percentage lowest meritValues

    if df_Theo.iloc[0]['rot'] == df_Theo.iloc[1]['rot'] == df_Theo.iloc[2]['rot'] == df_Theo.iloc[-1]['rot']: # rotation is fixed, don't plot it
        df_EE = df_Theo_EE.drop(columns=['rot', 'rot_err', 'logTeff', 'logL', 'logg'], errors='ignore') # make new dataframe without the spectroscopic info
        df = df_Theo.drop(columns=['rot', 'rot_err', 'logTeff', 'logL', 'logg'], errors='ignore') # make new dataframe without the spectroscopic info
        # Remove models in the error ellips from the regular dataframe.
        df = pd.merge(df,df_EE, indicator=True, how='outer', on=['Z', 'M', 'logD', 'fov', 'aov', 'Xc'], suffixes=[None, '_remove']).query('_merge=="left_only"').drop(['meritValue_remove', '_merge'], axis=1)

    else: # rotation was varied, include it in the plots
        df_EE = df_Theo_EE.drop(columns=['rot_err', 'logTeff', 'logL', 'logg'], errors='ignore') # make new dataframe without the spectroscopic info
        df = df_Theo.drop(columns=['rot_err', 'logTeff', 'logL', 'logg'], errors='ignore') # make new dataframe without the spectroscopic info
        # Remove models in the error ellips from the regular dataframe.
        df = pd.merge(df,df_EE, indicator=True, how='outer', on=['Z', 'M', 'logD', 'fov', 'aov', 'Xc'], suffixes=[None, '_remove']).query('_merge=="left_only"').drop(['meritValue_remove', 'rot_remove', '_merge'], axis=1)

    ax_dict={}  # dictionary of dictionaries, holding the subplots of the figure, keys indicate position (row, column) of the subplot
    nr_params = len(df.columns)-1
    for i in range(nr_params):
        ax_dict.update({i:{}})

    fig=plt.figure(figsize=(10,8))
    gs=GridSpec(nr_params,nr_params) # multiple rows and columns
    # proper label format on figures
    axis_labels_dict = {'rot': r'$\Omega_{\mathrm{rot}}$ [d$^{-1}$]' ,'M': r'M$_{\rm ini}$', 'Z': r'Z$_{\rm ini}$', 'logD':r'log(D$_{\rm env}$)', 'aov':r'$\alpha_{\rm CBM}$','fov':r'f$_{\rm CBM}$','Xc':r'$\rm X_c$'}

    if mark_best_model: min_index = df_EE['meritValue'].idxmin(axis='index', skipna=True)    # get the best model according to the point estimator

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
                    ax.hist( df_EE.iloc[:,nr_params-ix], bins=25, range=domain, density=False, cumulative=False, histtype='step' )

                else:
                    if len(values) > 1:
                        bin_half_width = (values[0]+values[1])/2-values[0]
                    else:
                        bin_half_width = 1E-3
                    bin_edges = [values[0]-bin_half_width]
                    for i in range(len(values)-1):
                        bin_edges.extend([(values[i]+values[i+1])/2])
                    bin_edges.extend([values[-1]+bin_half_width])
                    ax.hist( df_EE.iloc[:,nr_params-ix], bins=bin_edges, density=False, cumulative=False, histtype='step' )

                ax.tick_params(axis='y',left=False)
                continue

            im = ax.scatter(df.iloc[:,nr_params-ix], df.iloc[:,iy+1], c=np.log10(df.iloc[:,0]), cmap='Greys_r')
            im = ax.scatter(df_EE.iloc[:,nr_params-ix], df_EE.iloc[:,iy+1], c=np.log10(df_Theo_EE['meritValue']), cmap=CustomCMap)
            if mark_best_model: ax.scatter(df_EE.loc[min_index][nr_params-ix], df.loc[min_index][iy+1], color='white', marker = 'x')
            # Adjust x an y limits of subplots
            limit_adjust = (max(df.iloc[:,iy+1]) - min(df.iloc[:,iy+1]))*0.08
            ax.set_ylim( min(df.iloc[:,iy+1])-limit_adjust,  max(df.iloc[:,iy+1])+limit_adjust  )
            limit_adjust = (max(df.iloc[:,nr_params-ix]) - min(df.iloc[:,nr_params-ix])) *0.08
            ax.set_xlim( min(df.iloc[:,nr_params-ix])-limit_adjust, max(df.iloc[:,nr_params-ix])+limit_adjust )

    fig.align_labels()
    # add subplot in top right for Kiel or HRD
    ax_hrd = fig.add_axes([0.508, 0.65, 0.33, 0.33]) # X, Y, widht, height

    ax_hrd.set_xlabel(r'log(T$_{\mathrm{eff}}$ [K])', size=label_size)
    ax_hrd.tick_params(labelsize=label_size-4)
    ax_hrd.invert_xaxis()
    if logg_or_logL=='logg': ax_hrd.invert_yaxis()

    im = ax_hrd.scatter(df_Theo['logTeff'], df_Theo[logg_or_logL], c=np.log10(df_Theo['meritValue']), cmap='Greys_r')
    im_EE = ax_hrd.scatter(df_Theo_EE['logTeff'], df_Theo_EE[logg_or_logL], c=np.log10(df_Theo_EE['meritValue']), cmap=CustomCMap)
    ax_hrd.set_ylabel(f'{logg_or_logL[:-1]} {logg_or_logL[-1]}')
    if logg_or_logL == 'logL':
        ax_hrd.set_ylabel(r'log(L/L$_{\odot}$)', size=label_size)
    elif logg_or_logL == 'logg':
        ax_hrd.set_ylabel(r'$\log\,g$ [dex]', size=label_size)

    ax_hrd.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax_hrd.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax_hrd.tick_params(which='major', length=6)
    ax_hrd.tick_params(which='minor', length=4)

    # observations
    if n_sigma_spectrobox != None:
        Obs_dFrame  = pd.read_table(observations_file, delim_whitespace=True, header=0)
        # Observed spectroscopic error bar
        # To add the 1 and n-sigma spectro error boxes, calculate their width (so 2 and 2*n sigmas wide)
        width_logTeff_sigma= np.log10(Obs_dFrame['Teff'][0]+Obs_dFrame['Teff_err'][0]) - np.log10(Obs_dFrame['Teff'][0]-Obs_dFrame['Teff_err'][0])
        width_logTeff_nsigma= np.log10(Obs_dFrame['Teff'][0]+n_sigma_spectrobox*Obs_dFrame['Teff_err'][0]) - np.log10(Obs_dFrame['Teff'][0]-n_sigma_spectrobox*Obs_dFrame['Teff_err'][0])
        errorbox_1s = patches.Rectangle((np.log10(Obs_dFrame['Teff'][0]-Obs_dFrame['Teff_err'][0]),Obs_dFrame[logg_or_logL][0]-Obs_dFrame[f'{logg_or_logL}_err'][0]),
                    width_logTeff_sigma, 2*Obs_dFrame[f'{logg_or_logL}_err'][0],linewidth=1.7,edgecolor='cyan',facecolor='none')
        errorbox_ns = patches.Rectangle((np.log10(Obs_dFrame['Teff'][0]-n_sigma_spectrobox*Obs_dFrame['Teff_err'][0]), Obs_dFrame[logg_or_logL][0]-n_sigma_spectrobox*Obs_dFrame[f'{logg_or_logL}_err'][0]),
                    width_logTeff_nsigma, 2*n_sigma_spectrobox*Obs_dFrame[f'{logg_or_logL}_err'][0],linewidth=1.7,edgecolor='cyan',facecolor='none')
        ax_hrd.add_patch(errorbox_1s)
        ax_hrd.add_patch(errorbox_ns)
    if mark_best_model: ax_hrd.scatter(df_Theo_EE['logTeff'][min_index], df_Theo_EE[logg_or_logL][min_index], marker='x', color='white')

    # Add color bar
    cax = fig.add_axes([0.856, 0.565, 0.04, 0.415]) # X, Y, widht, height
    cbar= fig.colorbar(im, cax=cax, orientation='vertical')
    cax2 = fig.add_axes([0.856, 0.137, 0.04, 0.415]) # X, Y, widht, height
    cbar2= fig.colorbar(im_EE, cax=cax2, orientation='vertical', )

    if df_Theo_EE.shape[0]==1: # To prevent messing up colors due to automatic rescaling of colorbar
        im_EE.set_clim(np.log10(df_Theo_EE['meritValue']), np.log10(df_Theo_EE['meritValue'])*1.1)

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
    Path(fig_outputDir).mkdir(parents=True, exist_ok=True)
    fig.savefig(f'{fig_outputDir}{fig_title}.png', dpi=400)
    plt.clf()
    plt.close('all')

################################################################################

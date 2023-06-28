"""Extract frequencies from a grid of GYRE output, and generate period spacing series."""
import numpy as np
import pandas as pd
from pathlib import Path
import glob, logging
import multiprocessing
from functools import partial
from foam import support_functions as sf

logger = logging.getLogger('logger.ffg')

################################################################################
def generate_spacing_series(periods, errors=None):
    """
    Generate the period spacing series (delta P = p_(n+1) - p_n )
    ------- Parameters -------
    periods, errors (optional): list of floats
        Periods and their errors in units of days
    ------- Returns -------
    observed_spacings, observed_spacings_errors: tuple of lists of floats
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
def extract_frequency_grid(gyre_files, output_file='pulsationGrid.hdf', parameters=['rot', 'Z', 'M', 'logD', 'aov', 'fov', 'Xc'], nr_cpu=None):
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
    nr_cpu: int
        Number of worker processes to use in multiprocessing. The default 'None' will use the number returned by os.cpu_count().
    """
    MP_list = multiprocessing.Manager().list()    # make empty MultiProcessing listProxy

    # Glob all the files, then iteratively send them to a pool of processors
    summary_files = glob.iglob(gyre_files)
    with multiprocessing.Pool(nr_cpu) as p:
        extract_func = partial(all_freqs_from_summary, parameters=parameters)
        dictionaries = p.imap(extract_func, summary_files)
        for new_row in dictionaries:
            MP_list.append(new_row)   # Fill the listProxy with dictionaries for each read file

        df = pd.DataFrame(data=list(MP_list))
    # Sort the columns with frequencies by their radial order
    column_list = list(df.columns[:len(parameters)])
    column_list.extend(sorted(df.columns[len(parameters):]))
    df = df.reindex(column_list, axis=1)

    # Generate the directory for the output file and write the file afterwards
    Path(Path(output_file).parent).mkdir(parents=True, exist_ok=True)
    df.to_hdf(f'{output_file}', 'pulsgrid', format='table', mode='w')

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

    attributes, data = sf.read_hdf5(GYRE_summary_file)
    param_dict = sf.get_param_from_filename(GYRE_summary_file, parameters, values_as_float=True)

    for j in range(len(data['freq'])-1, -1, -1):    # Arrange increasing in radial order
        n_pg = data["n_pg"][j]
        if abs(n_pg) < 10:
            n_pg = f'{sf.sign(n_pg)}00{abs(n_pg)}'
        elif abs(n_pg) < 100:
            n_pg = f'{sf.sign(n_pg)}0{abs(n_pg)}'
        param_dict.update({f'n_pg{n_pg}':data['freq'][j][0]})

    return param_dict

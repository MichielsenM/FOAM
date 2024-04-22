"""Extract frequencies from a grid of GYRE output, and generate period spacing series."""

import glob
import logging
import multiprocessing
from functools import partial
from pathlib import Path

import pandas as pd

from foam import support_functions as sf

logger = logging.getLogger("logger.ffg")


################################################################################
def extract_frequency_grid(
    gyre_files, output_file="pulsationGrid.hdf", parameters=["rot", "Z", "M", "logD", "aov", "fov", "Xc"], nr_cpu=None
):
    """
    Extract frequencies from each globbed GYRE file and write them to 1 large file.

    Parameters
    ----------
    gyre_files: string
        String to glob to find all the relevant GYRE summary files.
    output_file: string
        Name (can include a path) for the file containing all the pulsation frequencies of the grid.
    parameters: list of strings
        List of parameters varied in the computed grid, so these are taken from the
        name of the summary files, and included in the 1 file containing all the info of the whole grid.
    nr_cpu: int
        Number of worker processes to use in multiprocessing.
        The default 'None' will use the number returned by os.cpu_count().
    """
    # make empty MultiProcessing listProxy
    mp_list = multiprocessing.Manager().list()

    # Glob all the files, then iteratively send them to a pool of processors
    summary_files = glob.iglob(gyre_files)
    with multiprocessing.Pool(nr_cpu) as p:
        extract_func = partial(all_freqs_from_summary, parameters=parameters)
        dictionaries = p.imap(extract_func, summary_files)
        for new_row in dictionaries:
            # Fill the listProxy with dictionaries for each read file
            mp_list.append(new_row)

        df = pd.DataFrame(data=list(mp_list))
    # Sort the columns with frequencies by their radial order
    column_list = list(df.columns[: len(parameters)])
    column_list.extend(sorted(df.columns[len(parameters) :]))
    df = df.reindex(column_list, axis=1)

    # Generate the directory for the output file and write the file afterwards
    Path(Path(output_file).parent).mkdir(parents=True, exist_ok=True)
    df.to_hdf(path_or_buf=output_file, key="pulsation_grid", format="table", mode="w")


################################################################################
def all_freqs_from_summary(gyre_summary_file, parameters):
    """
    Extract model parameters and pulsation frequencies from a GYRE summary file

    Parameters
    ----------
    gyre_summary_file: string
        path to the GYRE summary file
    parameters: list of strings
        List of input parameters varied in the computed grid,
        so these are read from the filename and included in returned line.

    Returns
    ----------
    param_dict: dict
        Dictionary containing all the model parameters and pulsation frequencies of the GYRE summary file.
    """

    _, data = sf.read_hdf5(gyre_summary_file)
    param_dict = sf.get_param_from_filename(gyre_summary_file, parameters, values_as_float=True)

    # Arrange increasing in radial order
    for j in range(len(data["freq"]) - 1, -1, -1):
        n_pg = data["n_pg"][j]
        if abs(n_pg) < 10:
            n_pg = f"{sf.sign(n_pg)}00{abs(n_pg)}"
        elif abs(n_pg) < 100:
            n_pg = f"{sf.sign(n_pg)}0{abs(n_pg)}"
        param_dict.update({f"n_pg{n_pg}": data["freq"][j][0]})

    return param_dict

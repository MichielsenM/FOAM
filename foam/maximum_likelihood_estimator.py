""" Functions to perform different kinds of maximum likelihood estimation for the models in a grid, and make correlation plots.
Note: The file with observations needs to hold temperature as Teff, although the analysis is done using the logTeff values."""

import logging
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from foam import build_optimised_pattern as bop
from foam import support_functions as sf

# Make a child logger of "logger" made in the top level script
logger = logging.getLogger("logger.mle_estimator")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def calculate_likelihood(
    theory_file,
    observables=None,
    merit_function=None,
    obs_path=None,
    star_name=None,
    fixed_params=None,
    grid_parameters=None,
):
    """
    Perform a maximum likelihood estimation using the provided type of merit function on the list of observables.
    Writes a data file with the values of the merit function and input parameters of each model.
    Can select and continue the analysis of nested grids through the keyword 'fixed_params'.

    Parameters
    ----------
    obs_path: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    theory_file: string
        Path to the hdf5 file with the theoretical model input parameters (first set of columns), frequency or period values (last set of columns),
        and possibly extra columns with additional observables (these columns should be in between the input parameters and frequency columns).
    observables: list of strings
        Which observables are included in the merit function.
        Must contain either 'f' (frequency), 'P' (period), or 'dP' (period-spacing) which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.
    merit_function: string
        The type of merit function to use. Currently supports "CS and "MD" ("chi-squared" and "mahalanobis distance").
    star_name: string
        Name of the star, used in file naming.
    fixed_params: dict
        Only select and analyse the part of the theoretical grid with the specified parameter values.
        The keys specify for which parameters only the specified value should be selected.
    grid_parameters: list of string
        List of the parameters in the theoretical grid.
    """
    if "f" in observables:
        observed_quantity = "frequency"
    elif "P" or "dP" in observables:
        observed_quantity = "period"
    # Read in the observed data and make an array of the observed observables
    obs_dataframe = pd.read_table(obs_path, delim_whitespace=True, header=0, index_col="index")
    obs, obs_err, file_suffix_observables = create_obs_observables_array(obs_dataframe, observables)

    # set the name of the output file
    _, tail = sf.split_line(Path(theory_file).stem, star_name)
    filename = f"{star_name}{tail}_{merit_function}_{file_suffix_observables}"

    # Theoretical grid data
    theory_dataframe = sf.get_subgrid_dataframe(theory_file, fixed_params)

    # get the interruptions in the pattern, absolute index in dataframe
    missing_absolute = np.where(theory_dataframe.columns.to_series().str.contains(f"{observed_quantity}_missing"))[0]
    # get the interruptions in the pattern, index relative within pulsations
    missing_relative = np.where(
        theory_dataframe.filter(like=f"{observed_quantity}")
        .columns.to_series()
        .str.contains(f"{observed_quantity}_missing")
    )[0]
    # Remove columns of missing frequencies
    theory_dataframe = theory_dataframe.drop(columns=theory_dataframe.columns[missing_absolute])
    # Adjust indices for removed lines of missing frequencies
    missing_indices = [missing_relative[i] - i for i in range(len(missing_relative))]

    thetas = np.asarray(theory_dataframe.filter(["rot"] + ["rot_err"] + grid_parameters))
    theory_puls = np.asarray(theory_dataframe.filter(like=f"{observed_quantity}"))

    # Make new list of theoretical models without entries with value -1
    new_theory = []
    new_thetas = []
    for i in range(len(theory_puls)):
        # ignore models where one of the freqs is -1
        if (theory_puls[i][0] != -1) and (theory_puls[i][-1] != -1):
            # make an array of the theoretical observables for each model
            theory_observables = create_theory_observables_array(theory_dataframe, i, observables, missing_indices)

            new_theory.append(theory_observables)
            new_thetas.append(thetas[i])
    neg_value = np.unique(np.where(theory_puls == -1)[0])
    logger.debug(f"File: {theory_file}")
    logger.debug(f"observables      : {observables}")
    logger.debug(f"ignored -1 freqs : {len(neg_value)}")
    logger.debug(f"original #models : {len(theory_puls)}")
    logger.debug(f"remaining #models: {len(new_theory)}")
    if len(neg_value) > 0:
        logger.warning(
            f"""{len(neg_value)} models were discarded due to mismatches during the selection of theoretical frequencies.
                        This is likely due to the frequency range used for the theoretical calculations being too narrow."""
        )
    new_theory = np.asarray(new_theory)
    new_thetas = np.asarray(new_thetas)

    # Dictionary containing different merit functions
    switcher = {"CS": merit_chi2, "MD": merit_mahalanobis}

    # get the desired function from the dictionary. Returns the lambda function if option is not in the dictionary.
    selected_merit_function = switcher.get(
        merit_function,
        lambda x, y, z: sys.exit(logger.error("invalid type of maximum likelihood estimator")),
    )
    merit_values = selected_merit_function(obs, obs_err, new_theory, fig_title=f"{filename}", star_name=star_name)

    # Combine values and save the results
    # add an additional column for MLE 'meritValues'
    combined_data = np.concatenate((np.matrix(merit_values).T, new_thetas), axis=1)
    # put the data in a pandas DataFrame
    df = pd.DataFrame(data=combined_data, columns=["meritValue"] + ["rot"] + ["rot_err"] + grid_parameters)
    df = pd.merge(
        df,
        theory_dataframe.drop(list(theory_dataframe.filter(regex=f"{observed_quantity}").columns), axis=1),
        how="inner",
        on=["rot", "rot_err"] + grid_parameters,
    )
    df.to_hdf(f"{os.getcwd()}/meritvalues/{filename}.hdf", "merit_values", format="table", mode="w")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def create_theory_observables_array(theory_dataframe, index, observables_in, missing_indices):
    """
    Create an array of theoretical observables.

    Parameters
    ----------
    theory_dataframe: pandas dataFrame
        DataFrame containing the theoretical periods or frequencies (as the last columns), along with any additional
        columns containing extra observables.
    index: int
        Row index in the dataFrame of the theoretical model to make the array for.
    observables_in: list of strings
        Which observables are included in the returned array.
        Must contain either 'f' (frequency), 'P' (period), or 'dP' (period-spacing) which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.
    missing_indices: list of int
        Contains the indices of the missing pulsations so that the period spacing pattern can be split around them.

    Returns
    ----------
    observables_out: numpy array, dtype=float
        The values of the specified observables for the model.
    """
    # Make a copy to leave the array handed to this function unaltered.
    observables = list(observables_in)
    if "P" in observables:
        # Cast to list before using asarray to prevent memory leak
        # add the periods to the output list
        observables_out = np.asarray(list(theory_dataframe.filter(like="period").loc[index]))
        observables.remove("P")

    elif "f" in observables:
        # Cast to list before using asarray to prevent memory leak
        # add the frequencies to the output list
        observables_out = np.asarray(list(theory_dataframe.filter(like="frequency").loc[index]))
        observables.remove("f")

    elif "dP" in observables:
        periods = np.asarray(theory_dataframe.filter(like="period").loc[index])
        observables_out = []
        for periods_part in np.split(periods, missing_indices):
            spacing, _ = bop.generate_spacing_series(periods_part)
            # switch back from seconds to days (so both P and dP are in days)
            observables_out = np.append(observables_out, np.asarray(spacing) / 86400)
        observables.remove("dP")

    # Add all other observables in the list from the dataFrame
    for observable in observables:
        observables_out = np.append(observables_out, np.asarray(theory_dataframe.loc[index, observable]))

    return observables_out


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def create_obs_observables_array(obs_dataframe, observables):
    """
    Create an array of the observed observables.

    Parameters
    ----------
    obs_dataframe: pandas dataFrame
        DataFrame containing the theoretical frequencies, periods, and any additional observables as columns, as well as columns with their errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    observables: list of strings
        Which observables are included in the returned array.
        Must contain either 'f' (frequency), 'P' (period), or 'dP' (period-spacing) which will be computed for the period pattern.
        Can contain any additional observables that are added as columns in both the file with observations and the file with theoretical models.

    Returns
    ----------
    observables_out: numpy array, dtype=float
        The values of the specified observables.
    observables_err_out: numpy array, dtype=float
        The errors on the values in observables_out.
    filename_suffix: string
        suffix for the filename, containing all the included observables, separated by '-'
    """
    # get the interruptions in the pattern
    missing_indices = np.where(obs_dataframe.index.isin(["f_missing"]))[0]
    # Adjust indices for removed lines of missing frequencies
    missing_indices = [missing_indices[i] - i for i in range(len(missing_indices))]
    # remove lines indicating missing frequencies (if they are present)
    if len(missing_indices) != 0:
        obs_dataframe = obs_dataframe.drop(index="f_missing")

    observables = list(observables)  # make a copy of the list, to not alter the one that was given to the function
    filename_suffix = ""

    if "P" in observables:
        observables_out = np.asarray(obs_dataframe["period"])
        observables_err_out = np.asarray(obs_dataframe["period_err"])
        filename_suffix = "P"
        observables.remove("P")

    elif "f" in observables:
        observables_out = np.asarray(obs_dataframe["frequency"])
        observables_err_out = np.asarray(obs_dataframe["frequency_err"])
        filename_suffix = "f"
        observables.remove("f")

    elif "dP" in observables:
        observables_out = []
        observables_err_out = []
        period = np.asarray(obs_dataframe["period"])
        period_err = np.asarray(obs_dataframe["period_err"])
        periods_parts = np.split(period, missing_indices)
        periods_err_parts = np.split(period_err, missing_indices)
        for periods, periods_err in zip(periods_parts, periods_err_parts):
            spacing, spacing_errs = bop.generate_spacing_series(periods, periods_err)
            # switch back from seconds to days (so both P and dP are in days)
            observables_out = np.append(observables_out, np.asarray(spacing) / 86400)
            observables_err_out = np.append(observables_err_out, np.asarray(spacing_errs) / 86400)

        filename_suffix = "dP"
        observables.remove("dP")

    if len(observables) > 0:
        filename_suffix += "+extra"
    # Add all other observables in the list from the dataFrame
    for observable in observables:
        logTeff = False
        if observable == "logTeff":
            # To read it as Teff from the observations datafile
            observable = "Teff"
            logTeff = True
        # since these columns have less entries, the dataFrame has NaNs in the empty rows
        obs = np.asarray(obs_dataframe[observable])
        # remove all entries that are NaN in the numpy array
        obs = obs[~np.isnan(obs)]
        obs_err = np.asarray(obs_dataframe[f"{observable}_err"])
        obs_err = obs_err[~np.isnan(obs_err)]

        # Convert observed Teff and error to log values
        if logTeff:
            logT = np.log10(obs)
            logT_err = obs_err / (obs * np.log(10))
            observables_out = np.append(observables_out, logT)
            observables_err_out = np.append(observables_err_out, logT_err)
        else:
            observables_out = np.append(observables_out, obs)
            observables_err_out = np.append(observables_err_out, obs_err)

    return observables_out, observables_err_out, filename_suffix


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def merit_chi2(Yobs, obs_err, YTheo, fig_title=None, star_name=None):
    """
    Calculate chi squared values for the given theoretical patterns

    Parameters
    ----------
    Yobs: numpy array, dtype=float
        observed values (period or frequency)
    obs_err: numpy array, dtype=float
        Errors on Yobs
    YTheo: numpy ndarray (2D), dtype=float
        Array of all theoretical patterns to calculate the chi squared value for.
    fig_title: None
        Should not be used in this function, but is to make it analogous to merit_mahalanobis()
        and enable the use of the lambda function.
    star_name: None
        Should not be used in this function, but is to make it analogous to merit_mahalanobis()
        and enable the use of the lambda function.
    Returns
    ----------
    chi2: numpy array, dtype=float
        chi squared values for the given theoretical values
    """
    chi2 = np.array([np.sum(((one_YTheo - Yobs) / obs_err) ** 2) for one_YTheo in YTheo])
    return chi2


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def merit_mahalanobis(Yobs, obs_err, YTheo, generate_output=True, fig_title=None, star_name=None):
    """
    Calculate mahalanobis distance (MD) values for the given theoretical patterns.

    Parameters
    ----------
    Yobs: numpy array, dtype=float
        observed values (period or frequency)
    obs_err: numpy array, dtype=float
        Errors on Yobs
    YTheo: numpy array of arrays of floats
        Array of all theoretical patterns to calculate the MD value for.
    generate_output: boolean
        Flag to write output and plot the variance-covariance matrix
    fig_title: string
        The name of the figure to be created.
    star_name: string
        The name of the analysed star, for file naming purposes.

    Returns
    ----------
    MD: numpy array, dtype=float
        Mahalanobis distances for the given theoretical patterns.
    """
    # Convert to matrix format (np.matrix is not recommended, use array and newaxis instead)
    YobsMat = np.array(Yobs)[np.newaxis].T
    YTheoMat = np.array(YTheo)[np.newaxis].T

    # Calculate the average on the theoretical values (e.g. frequencies)
    # over the entire grid. Returns a vector of dimension N x 1
    Yav = YTheoMat.mean(1)

    # Calculate the variance-covariance matrix
    # number of grid points
    q = np.shape(YTheo)[0]
    # number of observed values
    N = len(Yobs)
    v_matrix = np.zeros((N, N))
    for i in range(q):
        difference = np.subtract(YTheoMat[:, i], Yav)
        v_matrix += np.matmul(difference, difference.T)
    v_matrix = v_matrix / float(q - 1)

    # Include observational errors in the variance-covariance matrix
    v_matrix = v_matrix + np.diag(obs_err**2.0)
    # check if positive definite and make figure
    check_matrix(v_matrix, generate_output=generate_output, fig_title=fig_title, star_name=star_name)
    # Calculate Mahalanobis distances
    MD = np.zeros(q)
    v_matrix_inv = np.linalg.inv(v_matrix)
    for i in range(q):
        diff = YTheoMat[:, i] - YobsMat
        MD[i] = np.matmul(np.matmul(diff.T, v_matrix_inv), diff)[0][0]

    return MD


################################################################################
def check_matrix(v_matrix, generate_output=True, fig_title="Vmatrix", star_name=None):
    """
    Check the if the the eigenvalues of the Variance-covariance matrix are all positive,
    since this means the matrix is positive definite. Compute its determinant and condition number,
    and write them to a tsv file. Create and save a figure of the variance-covariance matrix.

    Parameters
    ----------
    v_matrix: numpy ndarray (2D), dtype=float
        Variance-covariance matrix
    output: boolean
        Flag to write output and plot the variance-covariance matrix
    fig_title: string
        The name of the figure to be created.
    star_name: string
        The name of the analysed star, for file naming purposes.
    """
    # If all eigenvalues are >0, it is positive definite
    if np.all(np.linalg.eigvals(v_matrix) > 0) == False:
        logger.error("V matrix is possibly not positive definite (since eigenvalues are not all > 0)")
        sys.exit(1)
    logger.debug(f"max(v_matrix) = {np.max(v_matrix)}")
    # multiply the matrix by the exponent of this, otherwise the determinant can be too small for the numerics
    kk = 10

    if generate_output is True:
        file_path = Path(f"{os.getcwd()}/V_matrix/{star_name}_determinant_conditionNr.tsv")
        file_path.parent.mkdir(parents=True, exist_ok=True)
        if not file_path.is_file():
            with file_path.open("w") as file:
                file.write(f"method \t ln(det(V)) \t condition_number \n")
        with file_path.open("a") as file:
            file.write(
                f"{fig_title} \t {np.log(np.linalg.det(np.exp(kk)*v_matrix))-kk*v_matrix.shape[0]:.2f} \t {np.linalg.cond(v_matrix):.2f} \n "
            )

        # Do *10^4 to get rid of small values, and put this in the colorbar label
        im = plt.imshow(v_matrix * 10**4, aspect="auto", cmap="Reds")
        plt.ylabel(rf"obs {v_matrix.shape[0]}  $\leftarrow$ obs 1", size=14)
        plt.xlabel(rf"obs 1 $\rightarrow$ obs {v_matrix.shape[0]} ", size=14)

        cbar = plt.colorbar(im)
        cbar.ax.set_ylabel(r"[d$^{2} 10^{-4}$]", rotation=90, labelpad=15, size=14)
        cbar.ax.tick_params(labelsize=14)
        plt.title("Variance covariance matrix")
        plt.tick_params(labelsize=14)
        plt.tight_layout()
        plt.savefig(f"{os.getcwd()}/V_matrix/{fig_title}.png")
        plt.clf()
        plt.close("all")

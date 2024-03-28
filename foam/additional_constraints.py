""" Extra constraints based on e.g. surface properties of the star, or a binary companion. """

import logging
import sys
from functools import partial
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger("logger.ac")


################################################################################
def surface_constraint(
    merit_values_file,
    observations_file=None,
    nsigma=3,
    constraint_companion=None,
    isocloud_grid_summary=None,
    surface_grid_file=None,
    free_parameters=["Z", "M", "logD", "aov", "fov", "Xc"],
    evolution_parameter="Xc",
    evolution_step=-1e-2,
):
    """
    Enforce an n-sigma constraint on the models based on the surface observables.
    Save this as a file with prefix indicating how many sigma the error box was.

    Parameters
    ----------
    merit_values_file: string
        Path to the hdf5 files with the merit function values and the surface info of the models in the grid.
    observations_file: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    nsigma: int
        How many sigma you want to make the interval to accept models.
    constraint_companion: dict
        Information on the companion star. Set to None to model single stars,
        or provide this to include binary constraints using isochrone-clouds.
    isocloud_grid_summary: dict
        Nested dictionary, the keys at its two levels are metallicity and mass.
        Holds the surface properties of the grid for the isochrone-cloud modelling per combination of metallicity-mass.
    surface_grid_file: string
        File with the surface properties and ages of the model-grid.
    free_parameters: list of strings
        List of all the parameters varied in the model grid.
    evolution_parameter: string
        Name of the parameter that is used to track the evolutionary steps of the model.
    evolution_step: float
        Change in the evolutionary parameter from one step to the next (negative if quantity decreases, e.g. central hydrogen content Xc)
    """
    obs_dataframe = pd.read_table(observations_file, delim_whitespace=True, header=0, index_col="index")
    dataframe_theory = pd.read_hdf(merit_values_file)

    if "Teff" in obs_dataframe.columns:
        dataframe_theory = dataframe_theory[
            dataframe_theory.logTeff < np.log10(obs_dataframe["Teff"][0] + nsigma * obs_dataframe["Teff_err"][0])
        ]
        dataframe_theory = dataframe_theory[
            dataframe_theory.logTeff > np.log10(obs_dataframe["Teff"][0] - nsigma * obs_dataframe["Teff_err"][0])
        ]
    if "logg" in obs_dataframe.columns:
        dataframe_theory = dataframe_theory[
            dataframe_theory.logg < obs_dataframe["logg"][0] + nsigma * obs_dataframe["logg_err"][0]
        ]
        dataframe_theory = dataframe_theory[
            dataframe_theory.logg > obs_dataframe["logg"][0] - nsigma * obs_dataframe["logg_err"][0]
        ]
    if "logL" in obs_dataframe.columns:
        dataframe_theory = dataframe_theory[
            dataframe_theory.logL < obs_dataframe["logL"][0] + nsigma * obs_dataframe["logL_err"][0]
        ]
        dataframe_theory = dataframe_theory[
            dataframe_theory.logL > obs_dataframe["logL"][0] - nsigma * obs_dataframe["logL_err"][0]
        ]

    if constraint_companion is not None:
        if (isocloud_grid_summary is None) or (surface_grid_file is None):
            logger.error(
                "Please supply a summary file for the isocloud grid and a path to the file with the grid surface properties and ages."
            )
            sys.exit()

        surface_grid_dataframe = pd.read_hdf(surface_grid_file)

        func = partial(
            enforce_binary_constraints,
            constraint_companion=constraint_companion,
            isocloud_grid_summary=isocloud_grid_summary,
            nsigma=nsigma,
            surface_grid_dataframe=surface_grid_dataframe,
            free_parameters=free_parameters,
            evolution_parameter=evolution_parameter,
            evolution_step=evolution_step,
        )
        indices_to_drop = dataframe_theory.apply(func, axis=1)
        for index_to_drop in indices_to_drop:
            # Check if index_to_drop equals itself to filter out NaN values
            if (index_to_drop is not None) and (index_to_drop == index_to_drop):
                dataframe_theory.drop(index_to_drop, inplace=True)

    output_file = f"{nsigma}sigmaBox_{merit_values_file}"
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    dataframe_theory.to_hdf(f"{output_file}", "surface_constrained_models", format="table", mode="w")


################################################################################
def get_age(
    model, df, free_parameters=["Z", "M", "logD", "aov", "fov", "Xc"], evolution_parameter="Xc", evolution_step=-1e-2
):
    """
    Get the age of the models one step older and younger than the provided model.

    Parameters
    ----------
    model: pandas series
        Parameters of the model.
    df: pandas dataFrame
        Dataframe with the model parameters and age (and surface info) of the theoretical models.
    free_parameters: list of strings
        List of all the parameters varied in the model grid.
    evolution_parameter: string
        Name of the parameter that is used to track the evolutionary steps of the model.
    evolution_step: float
        Change in the evolutionary parameter from one step to the next (negative if quantity decreases, e.g. central hydrogen content Xc)

    Returns
    ----------
    min_age, max_age: tuple of integers
        Age of the model one step younger and older than the provided model,
        these are the minimum and maximum age to accept models in the isochrone-cloud.
    """

    # copy to prevent deletion in the list outside this function
    params = list(free_parameters)
    params.remove(evolution_parameter)
    for param in params:
        df = df.loc[np.isclose(getattr(df, param), getattr(model, param))]

    model_evolution_attr = getattr(model, evolution_parameter)
    grid_evolution_attr = getattr(df, evolution_parameter)

    if abs(model_evolution_attr - max(grid_evolution_attr)) < abs(0.5 * evolution_step):
        min_age = 0
        max_age = int((df.loc[np.isclose(grid_evolution_attr, model_evolution_attr + evolution_step)].age).iloc[0])
    elif abs(model_evolution_attr - min(grid_evolution_attr)) < abs(0.5 * evolution_step):
        min_age = int((df.loc[np.isclose(grid_evolution_attr, model_evolution_attr - evolution_step)].age).iloc[0])
        age = int((df.loc[np.isclose(grid_evolution_attr, model_evolution_attr)].age).iloc[0])
        max_age = age + age - min_age
    else:
        min_age = int((df.loc[np.isclose(grid_evolution_attr, model_evolution_attr - evolution_step)].age).iloc[0])
        max_age = int((df.loc[np.isclose(grid_evolution_attr, model_evolution_attr + evolution_step)].age).iloc[0])
    return min_age, max_age


################################################################################
def enforce_binary_constraints(
    df_theory_row,
    constraint_companion=None,
    isocloud_grid_summary=None,
    nsigma=3,
    surface_grid_dataframe=None,
    free_parameters=["Z", "M", "logD", "aov", "fov", "Xc"],
    evolution_parameter="Xc",
    evolution_step=-1e-2,
):
    """
    Enforce an n-sigma constraint on the models based on
    spectroscopic observations of the binary companion employing isochrone-clouds.
    Assumes the same metallicity 'Z' for both primary and secondary,
    masses 'M' compatible with observed mass ratio 'q', and ages similar within 1 grid step.

    Parameters
    ----------
    df_theory_row: tuple, made of (int, pandas series)
        tuple returned from pandas.iterrows(), first tuple entry is the row index of the pandas dataFrame
        second tuple entry is a pandas series, containing a row from the pandas dataFrame.
        (This row holds model parameters, the merit function value, and surface properties.)
    constraint_companion: dict
        Information on the companion star, including surface parameters, mass ratio (q), the errors,
        and a boolean indicating whether the primary or secondary star is assumed pulsating and hence being modelled.
    isocloud_grid_summary: dict
        Nested dictionary, the keys at its two levels are metallicity and mass.
        Holds the surface properties of the grid for the isochrone-cloud modelling per combination of metallicity-mass.
    nsigma: int
        How many sigma you want to make the interval to accept models.
    surface_grid_dataframe: pandas DataFrame
        DataFrame with the surface properties and ages of the model-grid.
    free_parameters: list of strings
        List of all the parameters varied in the model grid.
    evolution_parameter: string
        Name of the parameter that is used to track the evolutionary steps of the model.
    evolution_step: float
        Change in the evolutionary parameter from one step to the next (negative if quantity decreases, e.g. central hydrogen content Xc)

    Returns
    ----------
    index: int or None
        Index of the dataframe that needs to be removed if binary constraints do not allow the model to remain.
        Returns None if the binary constraints do not discard the model.
    """
    model = df_theory_row
    min_age, max_age = get_age(
        model,
        surface_grid_dataframe,
        free_parameters=free_parameters,
        evolution_parameter=evolution_parameter,
        evolution_step=evolution_step,
    )
    q = constraint_companion["q"]
    q_err = constraint_companion["q_err"]
    if constraint_companion["primary_pulsates"]:
        m2_min = round(model.M * (q - q_err), 1)
        m2_max = round(model.M * (q + q_err), 1)
    else:
        m2_min = round(model.M / (q + q_err), 1)
        m2_max = round(model.M / (q - q_err), 1)

    isocloud_dict = isocloud_grid_summary[f"{model.Z}"]
    for key_mass, df in zip(isocloud_dict.keys(), isocloud_dict.values()):
        # Convert from string to float for the comparisons
        key_mass = float(key_mass)
        if key_mass < m2_min or key_mass > m2_max:
            # Only keep models that fall within mass range
            continue
        else:
            df = df[(df.star_age < max_age) & (df.star_age > min_age)]

            if df.shape[0] == 0:
                continue

            # Check for all provided constraints if the track passes through the uncertainty region
            if constraint_companion["Teff"] is not None:
                df = df[
                    df.log_Teff < np.log10(constraint_companion["Teff"] + nsigma * constraint_companion["Teff_err"])
                ]
                df = df[
                    df.log_Teff > np.log10(constraint_companion["Teff"] - nsigma * constraint_companion["Teff_err"])
                ]
            if constraint_companion["logg"] is not None:
                df = df[df.log_g < constraint_companion["logg"] + nsigma * constraint_companion["logg_err"]]
                df = df[df.log_g > constraint_companion["logg"] - nsigma * constraint_companion["logg_err"]]
            if constraint_companion["logL"] is not None:
                df = df[df.log_L < constraint_companion["logL"] + nsigma * constraint_companion["logL_err"]]
                df = df[df.log_L > constraint_companion["logL"] - nsigma * constraint_companion["logL_err"]]
            if df.shape[0] > 0:
                # If some models fall within the constraints, return None to not remove the model.
                return None

    return model.name

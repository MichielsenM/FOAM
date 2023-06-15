""" Extra constraints based on e.g. surface properties of the star, or a binary companion. """
import sys
from functools import partial
import numpy as np
import pandas as pd
from pathlib import Path
################################################################################
def surface_constraint(merit_values_file, observations_file=None, nsigma=3, constraint_companion=None, isocloud_grid_summary=None, surfaceGrid_file=None):
    """
    Enforce an n-sigma constraint on the models based on the surface observables.
    Save this as a file with prefix indicating how many sigma the error box was.
    ------- Parameters -------
    merit_values_file: string
        Path to the hdf5 files with the merit funtion values and the surface info of the models in the grid.
    observations_file: string
        Path to the tsv file with observations, with a column for each observable and each set of errors.
        Column names specify the observable, and "_err" suffix denotes that it's the error.
    nsigma: int
        How many sigmas you want to make the interval to accept models.
    constraint_companion: dict
        Information on the companion star. Set to None to model single stars,
        or provide this to include binary constraints using isochrone-clouds.
    isocloud_grid_directory: string
        Directory holding the grid for the isochrone-cloud modelling.
    surfaceGrid_file: string
        File with the surface properties and ages of the model-grid.
    """
    Obs_dFrame = pd.read_table(observations_file, delim_whitespace=True, header=0)
    df_Theo = pd.read_hdf(merit_values_file)

    if 'Teff' in Obs_dFrame.columns:
        df_Theo = df_Theo[df_Theo.logTeff < np.log10(Obs_dFrame['Teff'][0]+nsigma*Obs_dFrame['Teff_err'][0])]
        df_Theo = df_Theo[df_Theo.logTeff > np.log10(Obs_dFrame['Teff'][0]-nsigma*Obs_dFrame['Teff_err'][0])]
    if 'logg' in Obs_dFrame.columns:
        df_Theo = df_Theo[df_Theo.logg < Obs_dFrame['logg'][0]+nsigma*Obs_dFrame['logg_err'][0]]
        df_Theo = df_Theo[df_Theo.logg > Obs_dFrame['logg'][0]-nsigma*Obs_dFrame['logg_err'][0]]
    if 'logL' in Obs_dFrame.columns:
        df_Theo = df_Theo[df_Theo.logL < Obs_dFrame['logL'][0]+nsigma*Obs_dFrame['logL_err'][0]]
        df_Theo = df_Theo[df_Theo.logL > Obs_dFrame['logL'][0]-nsigma*Obs_dFrame['logL_err'][0]]

    if constraint_companion is not None:
        if (isocloud_grid_summary is None) or (surfaceGrid_file is None):
            logger.error('Please supply a directory for the isocloud grid and a path to the file with the grid surface properties and ages.')
            sys.exit()

        surfaceGrid_dataFrame = pd.read_hdf(surfaceGrid_file)

        func = partial(enforce_binary_constraints, constraint_companion=constraint_companion, isocloud_grid_summary=isocloud_grid_summary, nsigma=nsigma, surfaceGrid_dataFrame=surfaceGrid_dataFrame)
        indices_to_drop = df_Theo.apply(func, axis=1)
        for index_to_drop in indices_to_drop:
            if (index_to_drop is not None) and (index_to_drop==index_to_drop):
                df_Theo.drop(index_to_drop, inplace=True)

    outputFile = f'{nsigma}sigmaBox_{merit_values_file}'
    Path(outputFile).parent.mkdir(parents=True, exist_ok=True)
    df_Theo.to_hdf(f'{outputFile}', 'surface_constrained_models', format='table', mode='w')

################################################################################
def get_age(model, df):
    """
    Get the age of the models one step older and younger than the provided model.
    ------- Parameters -------
        model: pandas series
            Parameters of the model.
        df: pandas dataframe
            Dataframe with the model parameters and age (and surface info) of the theoretical models.

    ------- Returns -------
    min_age, max_age: tuple of integers
        Age of the model one sep younger and older than the procided model,
        these are the minimum and maximum age to accept models in the isochrone-cloud.
    """
    unique_xc=pd.unique(df[ np.isclose(df.Z, model.Z) ].Xc)

    if abs(model.Xc-max(unique_xc))<1E-4:
        min_age = 0
        max_age = int((df.loc[np.isclose(df.Z, model.Z) & np.isclose(df.M, model.M) & np.isclose(df.logD, model.logD) & np.isclose(df.aov, model.aov) & np.isclose(df.fov, model.fov) & np.isclose(df.Xc, round(model.Xc-0.01, 2))].age).iloc[0])
    elif abs(model.Xc-min(unique_xc))<1E-4:
        min_age = int((df.loc[np.isclose(df.Z, model.Z) & np.isclose(df.M, model.M) & np.isclose(df.logD, model.logD) & np.isclose(df.aov, model.aov) & np.isclose(df.fov, model.fov) & np.isclose(df.Xc, round(model.Xc+0.01, 2))].age).iloc[0])
        age     = int((df.loc[np.isclose(df.Z, model.Z) & np.isclose(df.M, model.M) & np.isclose(df.logD, model.logD) & np.isclose(df.aov, model.aov) & np.isclose(df.fov, model.fov) & np.isclose(df.Xc, round(model.Xc, 2))].age).iloc[0])
        max_age = age+age-min_age
    else:
        min_age = int((df.loc[np.isclose(df.Z, model.Z) & np.isclose(df.M, model.M) & np.isclose(df.logD, model.logD) & np.isclose(df.aov, model.aov) & np.isclose(df.fov, model.fov) & np.isclose(df.Xc, round(model.Xc+0.01, 2))].age).iloc[0])
        max_age = int((df.loc[np.isclose(df.Z, model.Z) & np.isclose(df.M, model.M) & np.isclose(df.logD, model.logD) & np.isclose(df.aov, model.aov) & np.isclose(df.fov, model.fov) & np.isclose(df.Xc, round(model.Xc-0.01, 2))].age).iloc[0])
    return min_age, max_age

################################################################################
def enforce_binary_constraints(df_Theo_row, constraint_companion=None, isocloud_grid_summary=None, nsigma=3, surfaceGrid_dataFrame=None):
    """
    Enforce an n-sigma constraint on the models based on spectoscopic observations of the binary companion employing isochrone-clouds.
    ------- Parameters -------
    df_Theo_row: tuple, made of (int, pandas series)
        tuple retruned from pandas.iterrows(), first tuple entry is the row index of the pandas dataFrame
        second tuple entry is a pandas series, containing a row from the pandas dataFrame.
        (This row holds model parameters, the meritfunction value, and surface properties.)
    constraint_companion: dict
        Information on the companion star, including surface parameters, mass ratio (q), the errors,
        and a boolean indicating whether the primary or secondary star is assumed pulsating and hence being modelled.
    isocloud_grid_directory: string
        Directory holding the grid for the isochrone-cloud modelling.
    nsigma: int
        How many sigmas you want to make the interval to accept models.
    surfaceGrid_dataFrame: pandas DataFrame
        DataFrame with the surface properteis and ages of the model-grid.
    ------- Returns -------
    index: int or None
        Index of the dataframe that needs to be removed if binary constraints do not allow the model to remain.
        Returns None if the binary constraints do not discard the model.
    """
    model = df_Theo_row
    min_age, max_age = get_age(model, surfaceGrid_dataFrame)
    q= constraint_companion['q']
    q_err= constraint_companion['q_err']
    if constraint_companion['primary_pulsates']:
        M2_min = round(model.M*(q-q_err), 1)
        M2_max = round(model.M*(q+q_err), 1)
    else:
        M2_min = round(model.M/(q+q_err), 1)
        M2_max = round(model.M/(q-q_err), 1)

    isocloud_dict = isocloud_grid_summary[f'{model.Z}']
    for key_Mass, df in zip(isocloud_dict.keys(), isocloud_dict.values()):
        key_Mass = float(key_Mass) #Convert from string to float for the comparisons
        if key_Mass < M2_min or key_Mass > M2_max:
            continue    # Only keep models that fall within mass range
        else:
            df = df[(df.star_age < max_age) & (df.star_age > min_age)]

            if df.shape[0] == 0:
                continue

            # Check for all provided constraints if the track passes through the uncertainty region
            if constraint_companion['Teff'] is not None:
                df=df[df.log_Teff < np.log10(constraint_companion['Teff']+nsigma*constraint_companion['Teff_err'])]
                df=df[df.log_Teff > np.log10(constraint_companion['Teff']-nsigma*constraint_companion['Teff_err'])]
            if constraint_companion['logg'] is not None:
                df=df[df.log_g < constraint_companion['logg']+nsigma*constraint_companion['logg_err']]
                df=df[df.log_g > constraint_companion['logg']-nsigma*constraint_companion['logg_err']]
            if constraint_companion['logL'] is not None:
                df=df[df.log_L < constraint_companion['logL']+nsigma*constraint_companion['logL_err']]
                df=df[df.log_L > constraint_companion['logL']-nsigma*constraint_companion['logL_err']]
            if df.shape[0] > 0: #If some models fall within the constraints, return None to not remove the model.
                return None

    return model.name

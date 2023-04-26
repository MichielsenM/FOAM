# Foam

A custom-made python package for FOrward Asteroseismic Modelling.
Contains useful functions to make figures, MESA/GYRE grids, and help with data analysis.

## Installation
<details>
 <summary>Installation instructions (click to expand) </summary>

 
Git clone this repository, and install using poetry (https://python-poetry.org/docs/) with command `poetry install` in the folder with the `pyproject.toml` file. This will install the package with all its dependencies, using the dependency versions as specified in the `poetry.lock` file. (The package will be installed in editable mode, so it will link the package to the original location, meaning any changes to the original package will be reflected directly in your environment.)

If you do not wish to use poetry, you could install by running `pip install .`. Note that this will install it as a package in your python environment, but in non-editable mode.

</details>

## Contents

1. `foam`: The FOrward Asteroseismic Modelling python package
2. `LICENSE`: GNU general public license
3. `poetry.lock`: List of dependencies and their exact versions.
4. `pyproject.toml`: Installation configuration file.
5. `tests`: pytest for the package

This repository is an open-source package, GNU-licensed, and any improvements provided by the users are well accepted. See GNU License in LICENSE.

### Author
Developed by Mathias Michielsen
```
mathias.michielsen@kuleuven.be
Instituut voor Sterrenkunde
KU Leuven, Belgium
```

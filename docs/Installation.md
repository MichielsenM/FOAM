---
layout: default
title: Installation
---
# Installation

Please note that FOAM was written and tested on Linux. The installation procedure was tested on MacOS, but the functionality was not as extensively tested as on Linux. No tests for Windows have been performed.

## Prerequisites
FOAM is a python package which requires python 3.9 or newer. The package requires <a href="https://gyre.readthedocs.io/en/stable/" target="_blank"> GYRE</a> to be installed as well. Please follow their <a href="https://gyre.readthedocs.io/en/stable/ref-guide/installation.html" target="_blank"> installation instructions</a>.

FOAM depends on a number of python packages which are listed in the <a href="https://github.com/MichielsenM/FOAM/tree/master/pyproject.toml" target="_blank"> pyproject.toml</a> file. These dependencies will automatically be installed when following the installation instructions below.

## Installing FOAM
### Installation using poetry (editable mode)
Poetry uses <a href="https://python-poetry.org/docs/" target="_blank"> poetry</a> for dependency management. Follow their installation instructions to install poetry on your system as well. If you for some reason do not wish to use poetry, you can alternatively install the package in non-editable mode using pip (instructions further below).

Once python and poetry are installed, git clone the <a href="https://github.com/MichielsenM/FOAM" target="_blank"> FOAM</a> repository to a location of your choice your machine. You can install it in your currently active python virtual environment using poetry with command `poetry install` in the top-level folder of the project where the `pyproject.toml` and `poetry.lock` files are located. This will install the package with all its dependencies. The exact versions of the dependencies are specified in the `poetry.lock` file so that everyone will get the same versions of the dependencies. The package will be installed in editable mode, so it will link the package to the original location where you cloned the repository, meaning any changes to the original package will be reflected directly in your environment.

### Installation without poetry (non-editable mode)
Git clone the <a href="https://github.com/MichielsenM/FOAM" target="_blank"> FOAM</a> repository to a location of your choice your machine. You can install the package in your currently active python virtual environment by running command `pip install .` in the top-level folder of the project where the `pyproject.toml` and `poetry.lock` files are located. Note that this will install it as a package in your python environment, but in non-editable mode. To be able to perform the pytests, you will need to run the command `pip install pytest tomli` as well.


### Test if installation was successful
To check if the package was installed successfully, you can run the following line in your terminal to print the package version: `python -m foam version`.
Additionally you can run the pytests via command `python -m foam test` to see if there are no errors.
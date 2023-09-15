---
layout: default
title: Installation
---
# Installation

## Prerequisites
FoAM is a python package which requires python 3.9 or newer, additionally the package uses <a href="https://python-poetry.org/docs/" target="_blank"> poetry</a> for dependency management. Follow their installation instructions to install poetry on your system as well.

## Installing FoAM
Once python and poetry are installed, git clone the <a href="https://github.com/MichielsenM/FoAM" target="_blank"> FoAM</a> repository to a location of your choice your machine. You can install it in your currently active python virtual environment using poetry with command `poetry install` in the top-level folder of the project where the `pyproject.toml` and `poetry.lock` files are located. This will install the package with all its dependencies. The exact versions of the dependencies are specified in the `poetry.lock` file so that everyone will get the same versions of the dependencies. The package will be installed in editable mode, so it will link the package to the original location where you cloned the repository, meaning any changes to the original package will be reflected directly in your environment.

### Installing without poetry
Although not recommended, if you for some reason do not wish to use poetry, you could install by running `pip install .` instead. Note that this will install it as a package in your python environment, but in non-editable mode.


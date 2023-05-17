---
layout: default
title: Walkthrough
---
# Walkthrough

## Setting up the directory
  In your project folder, make a subfolder for the star you want to model.
  The outputs of the modelling pipeline will be placed in this subfolder.
  Copy file `pipeline.py` from `foam/pipeline/` to it, this is the template script that can be used to run the parts of the pipeline. (See [Pipeline modules](./Pipeline.md) for more information.)

  The pipeline will initially create an additional folder, `grid_summary`, where it will store summary files for the theoretical model grid(s) when it is executed. This way the same summary files can be used to model multiple observed stars.

  The final directory structure, after following all steps in this walkthrough, will look similar to this.

  <details>
  <summary> Example directory structure (click to expand) </summary>
  (You can choose different names for all the files and folders, only folder `grid_summary` will be automatically generated and has a fixed name.)
  <pre>
  project_folder
  │
  └───grid_summary
  │
  └───star1
  │   │   data_KIC000.tsv
  │   │   pipeline.py
  |
  └───star2
  │   │   data_KIC001.tsv
  │   │   pipeline.py
  </pre>
  </details>


## The observational data
 Add a file with the observations in a .tsv format (tab-separated values). Its columns are the observables and their errors (suffix '\_err' for the error corresponding to a column). Frequencies, periods, and their errors are required for the full functionality of the pipeline. Teff (effective temperature), log(g) (surface gravity), log(L/Lsun) (luminosity), and their errors can be added to provide additional constraints on the asteroseismic solutions. Note that effective temperature should be added as Teff, but will be processed as logTeff (to be consistent with luminosity 'L' and surfacge gravity 'g' which are used in log scale as well).

 Frequencies are preferred to be listed decreasing in value, so increasing in period. The code can handle the other way around if the theoretical grid is also extracted in that way. (If the frequencies are ordered decreasing in value, the highest frequency parts in an interrupted pattern will be matched first, which is preferred since those have a lower relative uncertainty). One or more missing frequencies in the observed pattern can be indicated by 'f_missing' as index and 0 values for the frequency, period, and their errors. This will let the code model it as an interrupted period spacing pattern.

An example of what the structure of such a file with the observational data looks like is given in the table below.

 |    | frequency | frequency_err | period | period_err | Teff | Teff_err | logg | logg_err | logL | logL_err |
 |----|:---------:|--------------:|:------:|:----------:|:----:|:--------:|:----:|:--------:|:----:|:--------:|
 | f1 | 1.11 | 4e-5 | 0.9009 | 3e-5 | 15200 | 150 | 3.8 | 0.1 | 2.21 | 0.04 |
 | f2 | 1.04 | 5e-5 | 0.9615 | 4e-5 | | | | | | |
 | f3 | 0.98 | 2e-5 | 1.0204 | 1e-5 | | | | | | |
 | f_missing | 0 | 0 | 0 | 0 | | | | | | |
 | f4 | 0.87 | 2e-5 | 1.1494 | 1e-5 | | | | | | |

## The theoretical model grid
The pipeline will first create summary files that contain all the relevant information from the theoretical model grid regarding the pulsations and the general properties of the models.

To do this, the grid of theoretical models should adhere to a certain structure.
Starting off, your folder with model grids can contain multiple grids with different sets of physics. In our example below, we have two grids, each in their own folder, with different temperature structures. We'll call them `Diffusive` and `Peclet`.
For each grid, we will have a MESA grid of stellar equilibrium models, stored in `MESA_out`, and a GYRE grid with all the pulsation frequencies of those models, stored in `GYRE_out`.

#### MESA grid of stellar equilibrium models
All MESA output files are divided in subfolders per initial mass-metallicity combination. The naming scheme gives Zini and Mini, followed by their respective values, divided by an underscore.
In each of these folders, the MESA history files and profiles are divided in folders `history` and `profiles`. The individual filenames are constructed to hold the names of the varied input parameters followed by their numerical values. The different input parameters are divided by underscores.
In the example below, there were 5 parameters varied per MESA stellar evolution track (Z, M, logD, aov, fov).

The history files, ending with suffix '.hist', will thus be constructed like this: `Z0.008_M3.00_logD0.00_aov0.000_fov0.000.hist`
The history files should at least contain the following columns: 'star_age', 'log_Teff', 'log_L', 'log_g'

The profile files will have the same parameters as the history files, with one additional parameter to indicate their evolutionary stage. In our example, we use the central hydrogen mass fraction (Xc) for this purpose. The profile names, ending with '.prof', will hence be constructed like this: `Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.10.prof`
The profiles should at least contain the following information in their header: 'age', 'logTeff', 'logL', 'logg'

Note that the choice of parameters, and the way they are named, can be changed. The only requirements are that the parameter name is followed by it's numerical value, the different parameters are divided by underscores, and the files have the correct suffix. The parameters and names of your choice should then later be provided to the configuration of the pipeline. (See [Pipeline configuration](./Configuration.md) for more information.)

#### GYRE grid with stellar pulsations
For each MESA profile, a GYRE summary file needs to be calculated containing at least the output columns 'freq' and 'n_pg', and written in HDF format. These GYRE summary files are stored in folders indicating the rotation frequency (in cycles per day) and the mode ID (k,m) that were used during the GYRE computations (e.g. `rot0.63_k0m1`).

Similar to the MESA grid, the GYRE grid is divided in subfolders per initial mass-metallicity combination. The naming scheme gives Zini and Mini, followed by their respective values, divided by an underscore.

The naming scheme of the GYRE summary files should be similar to the MESA profiles, with one additional parameter 'rot' indicating the rotation frequency used in the GYRE computations. This results in names like the following:
`rot0.6304_Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.10.HDF`

<details>
<summary> Example grid directory structure (click to expand) </summary>
<pre>
Model_grids   
│
└───Diffusive
│   │
│   └───MESA_out
|   |   |
│   |   └───Zini0.008_Mini3.00
│   |   │   |
│   |   │   └───history
│   |   │   |   |   Z0.008_M3.00_logD0.00_aov0.000_fov0.000.hist
│   |   │   |   |   Z0.008_M3.00_logD0.00_aov0.000_fov0.005.hist
│   |   │   |   |   ...
│   |   │   |
│   |   │   └───profiles
│   |   │   |   |   Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.10.prof
│   |   │   |   |   Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.11.prof
│   |   │   |   |   ...
│   |   │
│   |   └───Zini0.008_Mini3.10
│   |   │   |
│   |   │   └───history
│   |   │   |   |   Z0.008_M3.10_logD0.00_aov0.000_fov0.000.hist
│   |   │   |   |   Z0.008_M3.10_logD0.00_aov0.000_fov0.005.hist
│   |   │   |   |   ...
│   |   │   |
│   |   │   └───profiles
│   |   │   |   |   Z0.008_M3.10_logD0.00_aov0.000_fov0.000_Xc0.10.prof
│   |   │   |   |   Z0.008_M3.10_logD0.00_aov0.000_fov0.000_Xc0.11.prof
│   |   │   |   |   ...
│   |   │
│   |   │   ...
│   |   
│   └───GYRE_out
│   |   │
│   |   └─── rot0.6304_k0m1
|   |   |   |
│   |   |   └───Zini0.008_Mini3.00
|   |   |   |   rot0.6304_Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.10.HDF
|   |   |   |   rot0.6304_Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.11.HDF
|   |   |   |   ...
|   |   |   |
│   |   |   └───Zini0.008_Mini3.10
|   |   |   |   rot0.6304_Z0.008_M3.10_logD0.00_aov0.000_fov0.000_Xc0.10.HDF
|   |   |   |   rot0.6304_Z0.008_M3.10_logD0.00_aov0.000_fov0.000_Xc0.11.HDF
|   |   |   |   ...
|   |   |   |
│   |   |   └───...
│   |   │
│   |   └─── rot1.13_k1m0
|   |   |   |
│   |   |   └───Zini0.008_Mini3.00
|   |   |   |   rot1.13_Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.10.HDF
|   |   |   |   rot1.13_Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.11.HDF
|   |   |   |   ...
|   |   |   |
│   |   |   └───Zini0.008_Mini3.10
|   |   |   |   rot1.13_Z0.008_M3.10_logD0.00_aov0.000_fov0.000_Xc0.10.HDF
|   |   |   |   rot1.13_Z0.008_M3.10_logD0.00_aov0.000_fov0.000_Xc0.11.HDF
|   |   |   |   ...
|   |   |   |
│   |   |   └───...
│
└───Peclet
│   │
│   └───MESA_out
│   |   │
│   |   └─── ...
│   │
│   └───GYRE_out
│   |   │
│   |   └───...

</pre>
</details>

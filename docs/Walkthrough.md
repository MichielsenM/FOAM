---
layout: default
title: Walkthrough
---
# Walkthrough

## Setting up the directory
  In your project folder, make a subfolder for the star you want to model.
  The outputs of the modelling pipeline will be placed in this subfolder.
  Copy file <a href="https://github.com/MichielsenM/FOAM/blob/master/foam/pipeline/pipeline.py" target="_blank"> pipeline.py</a> from folder `foam/pipeline/` to it, this is the template script that can be used to run the parts of the pipeline. (See [Pipeline modules](./Pipeline.md) for more information.)

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
 Add a file with the observations in a .tsv format (tab-separated values). Its columns are the observables and their errors (suffix '\_err' for the error corresponding to a column). Frequencies, periods, and their errors are required for the full functionality of the pipeline. Additional observables such as Teff (effective temperature), log(g) (surface gravity), log(L/Lsun) (luminosity), and their errors can be added to provide additional constraints on the asteroseismic solutions. They can also be completely omitted from this file if desired. Note that effective temperature should be added as Teff, but will be processed as logTeff (to be consistent with luminosity 'L' and surfacge gravity 'g' which are used in log scale as well).
 Furthermore, any additional observables that you want to include in the merit function can be added in this file.

 Frequencies are preferred to be listed decreasing in value, so increasing in period. The pipeline can handle the other way around if the theoretical grid is also extracted in that way. (If the frequencies are ordered decreasing in value, the highest frequency parts in an interrupted pattern will be matched first, which is preferred since those have a lower relative uncertainty). One or more missing frequencies in the observed pattern can be indicated by 'f_missing' as index and 0 values for the frequency, period, and their errors. This will set the pipeline to model it as an interrupted period spacing pattern. (One entry of 'f_missing' denotes an interruption in the pattern due to one or more missing frequencies.) 
Keep in mind the occurence of an isolated frequency in between missing frequencies (f_(n-1), f_missing, f_n, f_missing, f_(n+1), f_(n+2)...) will not be taken into account in the merit function when period spacings are used as seismic observables (`observable_seismic=['dP']`, see [Pipeline configuration](./Configuration.md)), since they won't have a period spacing associated with them.

An example of what the structure of such a file with the observational data looks like is given in the table below.

  <details>
  <summary> Example observational data (click to expand) </summary>
  You can add additional observables that you want to include in the merit function to the observations by adding extra columns like the ones for logL and logL_err.
  (See also <a href="https://github.com/MichielsenM/FOAM/blob/master/example_setup/KIC7760680/data_KIC7760680.tsv" target="_blank"> data_KIC7760680.tsv </a> from the example setup explained further down on this page.)

  <table>
    <tr>
      <th>index</th>
      <th>frequency</th>
      <th>frequency_err</th>
      <th>period</th>
      <th>period_err</th>
      <th>Teff</th>
      <th>Teff_err</th>
      <th>logg</th>
      <th>logg_err</th>
      <th>logL</th>
      <th>logL_err</th>
    </tr>
    <tr>
      <td>f1</td>
      <td>1.11</td>
      <td>4e-5</td>
      <td>0.9009 </td>
      <td>3e-5 </td>
      <td>15200 </td>
      <td>200 </td>
      <td>3.8 </td>
      <td>0.1 </td>
      <td>2.21 </td>
      <td>0.04 </td>
    </tr>
    <tr>
      <td>f2 </td>
      <td>1.04 </td>
      <td>5e-5 </td>
      <td>0.9615 </td>
      <td>4e-5 </td>
    </tr>
    <tr>
      <td>f3 </td>
      <td>0.98 </td>
      <td>2e-5 </td>
      <td>1.0204 </td>
      <td>1e-5 </td>
    </tr>
    <tr>
      <td>f_missing </td>
      <td>0 </td>
      <td>0 </td>
      <td>0 </td>
      <td>0 </td>
    </tr>
    <tr>
      <td>f4 </td>
      <td>0.87 </td>
      <td>2e-5 </td>
      <td>1.1494 </td>
      <td>1e-5 </td>
    </tr>      
  </table> 
  </details>

## The theoretical model grid
The pipeline will first create summary files that contain all the relevant information from the theoretical model grid regarding the pulsations and the general properties of the models. 
To do this, the grid of theoretical models should adhere to a certain structure.
Starting off, your folder with model grids can contain multiple grids with different sets of physics. In our example below, we have two grids, each in their own folder, with different temperature structures. We'll call them `Diffusive` and `Peclet`.
For each grid, we will have a <a href="https://docs.mesastar.org/en/latest/index.html" target="_blank"> MESA</a> grid of stellar equilibrium models, stored in `MESA_out`, and a <a href="https://gyre.readthedocs.io/en/stable/" target="_blank"> GYRE</a> grid with all the pulsation frequencies of those models, stored in `GYRE_out`.

The current implementation is made for a grid of stellar equilibrium models computed by <a href="https://docs.mesastar.org/en/latest/index.html" target="_blank"> MESA</a>, whose pulsation frequencies are computed with <a href="https://gyre.readthedocs.io/en/stable/" target="_blank"> GYRE</a>. Some adjustments could be made to apply the modelling pipeline to grids computed by other codes. Depending on the structure of the output files from these other codes, the adjustments would likely include (but possibly not be limited to) the functions `read_mesa_file` and `info_from_profiles` from module `functions_for_mesa`, and `all_freqs_from_summary` from module `functions_for_gyre`.

#### MESA grid of stellar equilibrium models
All MESA output files are divided in subfolders per initial mass-metallicity combination. The naming scheme gives Zini and Mini, followed by their respective values, divided by an underscore.
In each of these folders, the MESA history files and profiles are divided in folders `history` and `profiles`. The individual filenames are constructed to hold the names of the varied input parameters followed by their numerical values. The different input parameters are divided by underscores.
In the example below, there were 5 parameters varied per MESA stellar evolution track (Z, M, logD, aov, fov).

The history files, ending with suffix '.hist', will thus be constructed like this: `Z0.008_M3.00_logD0.00_aov0.000_fov0.000.hist`.
The history files should at least contain the following columns: 'star_age', 'log_Teff', 'log_L', 'log_g'.

The profile files will have the same parameters as the history files, with one additional parameter to indicate their evolutionary stage. In our example, we use the central hydrogen mass fraction (Xc) for this purpose. The profile names, ending with '.prof', will hence be constructed like this: `Z0.008_M3.00_logD0.00_aov0.000_fov0.000_Xc0.10.prof`.
The profiles should at least contain the data column 'log_g', and the following information in their header: 'star_age', 'Teff', 'photosphere_L'.
Furthermore, any additional observables that you want to include in the merit function (and that hence were added to the file with observational data will) need to be present in the MESA profile header information as well.

Note that the choice of parameters, the way they are named, and how many parameters are included, can all be changed. The only requirements are that the parameter name is followed by it's numerical value, the different parameters are divided by underscores, and the files have the correct suffix ('.hist' or '.prof'). The parameters and names of your choice should then later be provided to the configuration of the pipeline. (See [Pipeline configuration](./Configuration.md) for more information.)

#### GYRE grid with stellar pulsations
For each MESA profile, a GYRE summary file needs to be calculated containing at least the output columns 'freq' and 'n_pg', and written in HDF format. These GYRE summary files are stored in folders indicating the rotation frequency (in cycles per day) and the mode ID (k,m) that were used during the GYRE computations (e.g. `rot0.63_k0m1`).

Similar to the MESA grid, the GYRE grid is divided in subfolders per initial mass-metallicity combination. The naming scheme gives Zini and Mini, followed by their respective values, divided by an underscore. (But just as for the MESA grid, you can choose the naming of these parameters yourself.)

The naming scheme of the GYRE summary files should be the same as for the MESA profiles, with one additional parameter 'rot' indicating the rotation frequency used in the GYRE computations. This results in names like the following:
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


## Example setup

A practical example of the setup can be found in the following directory <a href="https://github.com/MichielsenM/FOAM/tree/master/example_setup" target="_blank"> example_setup</a>.
The folder `KIC7760680` contains an example of a file with observational data (<a href="https://github.com/MichielsenM/FOAM/blob/master/example_setup/KIC7760680/data_KIC7760680.tsv" target="_blank"> data_KIC7760680.tsv </a> ), and two setups for the modelling pipeline. One for modelling the full grid, and one setup to perform the modelling for a subgrid of nested models (in this case for both Convective Boundary Mixing parameters fixed at 0). Some print statements have been included to keep you informed how far the pipeline has progressed, but these can easily be removed or adjusted.

If you calculated your own theoretical grids, you should set their path properly in the pipeline files with the function argument `grid_parent_directory`. However, for this walkthrough the folder `grid_summary` already contains the outputs of step 0 of the modelling pipeline. You can unpack these compressed files using terminal commands `tar -xf surfaceGrid_DO.tar.xz` and `tar -xf pulsationGrid_DO_rot0.4805_k0m1.tar.xz`.
The theoretical grid summary that is included here is the radiative grid from <a href="https://doi.org/10.1051/0004-6361/202039926" target="_blank"> Michielsen et al. (2021)</a>. If you were to run the pipelines, you would find similar results to the paper, but slight differences would be present. These differences occur since the near-core rotation rate of the star was assumed to be fixed in the paper. Meanwhile the modelling pipeline has undergone numerous updates, one of which is the optimisation of the stellar rotation rate for each theoretical equilibrium model separately. 
This optimisation is explained in more detail in <a href="https://arxiv.org/abs/2309.13123" target="_blank"> Michielsen et al. (2023)</a>.

Starting the modelling with this pre-made setup is now as easy as running `python pipeline.py` from your terminal (after making sure the right python environment is activated).
Note that you will see warnings `WARNING  file already existed: ../grid_summary/surfaceGrid_DO.hdf` and `WARNING  file already existed: ../grid_summary/pulsationGrid_DO_rot0.4805_k0m1.hdf`. This is normal since the code detects that the summary files are already present, and don't need to be extracted a second time. Step 0 of the modelling pipeline will therefore skip creating these files and notify you that they already existed. If you want to speed up running the example setup, you can add `pattern_methods = ['provided-pulsation', 'highest-frequency']` to the pipeline configuration options. This will then enable only 2 out of the 3 options for constructing the theoretical pulsation patterns. The method `chisq_longest_sequence` will be left out, since this one is far slower than the others and will add a considerable amount of time to step 1 of the pipeline. It will of course remove some of the time in the other steps as well. Alternatively, you can also increase `nr_cpu` depending on your computational setup.

As explained before, the expected results are slightly different from those published by <a href="https://doi.org/10.1051/0004-6361/202039926" target="_blank"> Michielsen et al. (2021)</a>. The expected results of this pipeline are therefore included in directory <a href="https://github.com/MichielsenM/FOAM/tree/master/example_setup/expected_outcome" target="_blank"> expected_outcome</a>. This folder contains the results for modelling the full grid and the nested subgrid. It contains the output of the last step (step 7), which is a table with the best model parameters for each combination of the chosen theoretical grid, seismic observables, and pattern construction methods. One table holds these for a chi-squared (CS) merit function, whilst the other holds these for a Mahalanobis Distance (MD) merit function. Furthermore some of the figures created by step 6 of the pipeline have been included as well. (See [Pipeline modules](./Pipeline.md) for more information on the output of the different steps in the modelling.)


### Example setup interrupted pattern
Another example is provided for the case of an interrupted pattern with some frequencies missing.
For this setup, we make a new list of frequencies, <a href="https://github.com/MichielsenM/FOAM/blob/master/example_setup/KIC7760680_missing_freq/data_KIC7760680_missing_freq.tsv" target="_blank"> data_KIC7760680_missing_freq.tsv </a>, where we presumed that frequencies f4, f6, f18, and f19 from the original list <a href="https://github.com/MichielsenM/FOAM/blob/master/example_setup/KIC7760680/data_KIC7760680.tsv" target="_blank"> data_KIC7760680.tsv </a> are now missing. Additionally, we removed the luminosity information from the data file. Some notable differences in the configuration of the pipeline are in the settings for `pattern_starting_pulsation`, `N_pattern_parts`, and the different number for `N_periods`. Although not necessary to run the pipeline, we changed `observable_seismic` from period spacings ('dP') to periods ('P') for demonstrative purposes. The expected results of this pipeline are included in directory <a href="https://github.com/MichielsenM/FOAM/tree/master/example_setup/expected_outcome_missing_freq" target="_blank"> expected_outcome_missing_freq</a>. The cornerplots that originally had an HR-diagram in the top right corner now show a Kiel diagram instead, since we removed the luminosity information from the observational data.
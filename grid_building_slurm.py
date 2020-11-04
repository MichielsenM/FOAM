"""Functions for building MESA and GYRE grids on the SLURM framework."""
# from PyPulse import grid_building_slurm as gbs
import numpy as np
import glob, os, sys
from pathlib import Path
from shutil import copyfile
import logging
from . import my_python_functions as mypy
from . import functions_for_gyre as ffg

logger = logging.getLogger('logger.gbs')
################################################################################
def write_bash_submit( jobname, list, bash_submit_list, walltime=1440, memory=3000, cpu=1,
                       bash_template = os.path.expandvars('$CONDA_PREFIX/lib/python3.7/site-packages/PyPulse/templates/SLURM_submit_list_template.sh')
                      ):
    """
    Write a bash script to sumbit a list of jobs to SLURM.
    ------- Parameters -------
    jobname, list, bash_submit_list: string
        name of the job, path to the list to submit, name of the file for the bash script
    bash_template: string
        path to the template of the bash script to read and modify
    walltime, memory, cpu: int
        specifying walltime (minutes), memory (mb) and number of cpus per task
    """
    with open(bash_template, 'r') as f:
            lines = f.readlines()
    replacements = {'JOBNAME' : f'{jobname}',
                    'WALLTIME': f'{walltime}',
                    'MEMORY'  : f'{memory}',
                    'CPU'     : f'{cpu}',
                    'LIST'    : f'{list}' }

    new_lines = []
    for line in lines:
        new_line = line
        for key in replacements:
                if (replacements[key] != ''):
                        new_line = new_line.replace(key, replacements[key])
        new_lines.append(new_line)
    logger.info(f'write {bash_submit_list}')

    with open(bash_submit_list, 'w') as f:
        f.writelines(new_lines)
    return

################################################################################
def make_mesa_setup(setup_directory=f'{os.getcwd()}/MESA_setup', output_dir=f'{os.getcwd()}/MESA_out',
                    work_dir=f'{os.getcwd()}/MESA_work_dir', Z_ini_list=[0.014], M_ini_list=[1], log_Dmix_list=[1], aov_list=[0], fov_list=[0]):
    """
    Construct a setup for a MESA grid with job lists to run on e.g. SLURM, and bash scripts to run each job list.
    ------- Parameters -------
    setup_directory, output_dir, work_dir: string
        paths to the directory where the bash setup is being made, to the directory where the MESA output will be stored, and to the MESA work directory.
    Z_ini_list, M_ini_list, log_Dmix_list, fov_list, aov_list: list of floats
        Lists of the parameter values to be computed in the grid.
    """
    if not os.path.exists(work_dir):
        logger.error(f'Specified MESA work directory does not exist: {work_dir}')
        sys.exit()
    Path(f'{setup_directory}/submit-lists-scripts').mkdir(parents=True, exist_ok=True)
    Path(f'{setup_directory}/lists').mkdir(parents=True, exist_ok=True)

    params_to_run = []
    lines_run_all_bash = []
    lines_run_all_bash.append('#!/bin/bash \n')

    index = 0
    parameter_combinations = len(Z_ini_list)*len(M_ini_list)*len(log_Dmix_list)*len(aov_list)*len(fov_list)

    for Z_ini in Z_ini_list:
        Z_ini = f'{Z_ini:.3f}'
        for M_ini in M_ini_list:
            M_ini = f'{M_ini:.2f}'
            for log_Dmix in log_Dmix_list:
                log_Dmix = f'{log_Dmix:.2f}'
                for a_ov in aov_list:
                    a_ov = f'{a_ov:.3f}'
                    for f_ov in fov_list:
                        f_ov = f'{f_ov:.3f}'

                        params_to_run.append(f'{setup_directory}/run_MESA.sh {Z_ini} {M_ini} {log_Dmix} {a_ov} {f_ov} {output_dir} {work_dir} \n')
                        index += 1
                        if index%1000 == 0 or index == parameter_combinations:
                            list_nr = int(np.floor(index/1000))+1
                            list_to_run_path = f'{setup_directory}/lists/list{list_nr}'
                            write_bash_submit(f'MESA{list_nr}', list_to_run_path, f'{setup_directory}/submit-lists-scripts/submit_list{list_nr}.sh')
                            # logger.info(f'write {list_to_run_path}')
                            with open(list_to_run_path, 'w') as f:
                                f.writelines(params_to_run)
                            params_to_run = []

                            if index%1000 == 0:
                                lines_run_all_bash.append(f'sbatch --array=1-1000 {setup_directory}/submit-lists-scripts/submit_list{list_nr}.sh \n')
                            else:
                                lines_run_all_bash.append(f'sbatch --array=1-{parameter_combinations%1000} {setup_directory}/submit-lists-scripts/submit_list{list_nr}.sh \n')


    lines_run_all_bash = sorted(lines_run_all_bash)     #A bash script to submit all other bash scripts
    with open(f'{setup_directory}/submit_all_bash.sh', 'w') as fobj:
        fobj.writelines(lines_run_all_bash)

    copyfile(os.path.expandvars('$CONDA_PREFIX/lib/python3.7/site-packages/PyPulse/templates/run_MESA.sh'), f'{setup_directory}/run_MESA.sh')
    return

################################################################################
def make_gyre_setup(setup_directory=f'{os.getcwd()}/GYRE_setup', output_dir=f'{os.getcwd()}/GYRE_out', mesa_dir=f'{os.getcwd()}/MESA_out',
                    npg_min=-50, npg_max=-1, azimuthal_order=1, degree=1, omega_rot=[0.0], unit_rot = 'CYC_PER_DAY', rotation_frame='INERTIAL'):
    """
    Construct a setup for a GYRE grid with job lists to run on e.g. SLURM, and bash scripts to run each job list.
    GYRE inlists and jobs will be created for each pulsation file found in the MESA directory.
    Scanning ranges in GYRE inlists will be set based on desired n_pg range and Asymptotic_dP in the MESA hist file.
    ------- Parameters -------
    setup_directory, output_dir, mesa_dir: string
        paths to the directory where the bash setup is being made, to the directory where the GYRE output will be stored, and to the MESA output directory.
    npg_min, npg_max, degree, azimuthal_order : int
        Quantum numbers of the modes to be calculated. Minimum and maximum radial order, degree (l) and azimuthal order (m)
    omega_rot: list of float
        rotation frequency of the model
    unit_rot: string
        unit of the rotation frequency, can be CYC_PER_DAY or CRITICAL (roche critical)
    rotation_frame: string
        rotational frame of reference for the pulsation freqencies
    """
    if not os.path.exists(mesa_dir):
        logger.error(f'Specified MESA output directory does not exist: {mesa_dir}')
        sys.exit()
    Path(f'{setup_directory}/submit-lists-scripts').mkdir(parents=True, exist_ok=True)
    Path(f'{setup_directory}/lists').mkdir(parents=True, exist_ok=True)
    Path(f'{setup_directory}/inlists').mkdir(parents=True, exist_ok=True)
    Path(f'{output_dir}').mkdir(parents=True, exist_ok=True)

    lines_run_all_bash = []
    lines_run_all_bash.append('#!/bin/bash \n')
    gyre_files = glob.glob(mesa_dir + '/*/gyre/*.GYRE' )

    lines_to_run = []

    for rotation in omega_rot:
        for index, file_path in enumerate(gyre_files, start=1): # Start the index count from 1 instead of from 0, but still loop through all files
            path, filename = file_path.rsplit('/',1)
            param_dict = mypy.get_param_from_filename(file_path, ['M', 'Z', 'logD', 'aov', 'fov', 'Xc'])
            output_dir_Z = f'{output_dir}/rot{rotation}/Zini{param_dict["Z"]}/{filename[:filename.rfind(".")]}'

            f_min, f_max = ffg.calc_scanning_range(file_path, npg_min=npg_min, npg_max=npg_max, l=degree, m=azimuthal_order, omega_rot=rotation, unit_rot=unit_rot, rotation_frame=rotation_frame)
            inlist_to_write = f'{setup_directory}/inlists/rot{rotation}_{filename[:-5]}.in'
            write_gyre_inlist(inlist_to_write, file_path, npg_min=npg_min,npg_max=npg_max, freq_min=f_min, freq_max=f_max, omega_rot=rotation, unit_rot=unit_rot, rotation_frame=rotation_frame)
            lines_to_run.append(f'{setup_directory}/run_GYRE.sh {param_dict["Z"]} {param_dict["M"]} {param_dict["logD"]} {param_dict["aov"]} {param_dict["fov"]} {param_dict["Xc"]} {inlist_to_write} {output_dir_Z} \n')

            if index%1000 == 0 or index==len(gyre_files):
                list_nr = int(np.floor(index/1000))+1
                list_to_run_path = f'{setup_directory}/lists/rot{rotation}_list{list_nr}'
                write_bash_submit(f'GYRE{list_nr}', list_to_run_path, f'{setup_directory}/submit-lists-scripts/rot{rotation}_submit_list{list_nr}.sh', walltime=360, memory=500)
                with open(list_to_run_path, 'w') as f:
                    f.writelines(lines_to_run)
                lines_to_run = []

                if index%1000 == 0:
                    lines_run_all_bash.append(f'sbatch --array=1-1000 {setup_directory}/submit-lists-scripts/rot{rotation}_submit_list{list_nr}.sh \n')
                else:
                    lines_run_all_bash.append(f'sbatch --array=1-{index%1000} {setup_directory}/submit-lists-scripts/rot{rotation}_submit_list{list_nr}.sh \n')

    lines_run_all_bash = sorted(lines_run_all_bash)     #A bash script to submit all other bash scripts
    with open(f'{setup_directory}/submit_all_bash.sh', 'w') as fobj:
        fobj.writelines(lines_run_all_bash)

    copyfile(os.path.expandvars('$CONDA_PREFIX/lib/python3.7/site-packages/PyPulse/templates/run_GYRE.sh'), f'{setup_directory}/run_GYRE.sh')
    return
################################################################################
def write_gyre_inlist( gyre_in_file, mesa_pulsation_file, gyre_summary_file='',
                       freq_min=0.01, freq_max=10, rotation_frame='INERTIAL',
                       npg_min=-50,npg_max=-1, omega_rot=0.0, unit_rot = 'CYC_PER_DAY', azimuthal_order=1, degree=1,
                       gyre_base_file = os.path.expandvars('$CONDA_PREFIX/lib/python3.7/site-packages/PyPulse/templates/gyre_template.in')
                      ):
    """
    Write gyre inlists based upon a given template
    ------- Parameters -------
    gyre_in_file, mesa_pulsation_file, gyre_summary_file: strings
        paths to the GYRE inlist to be constructed, the MESA pulsation file to be read by the GYRE run, and GYRE output summary file
    freq_min_inertial, freq_max_inertial, omega_rot: float
        minimum and maximum frequency of the range to scan, rotation frequency in cyc/day
    rotation_frame: string
        gridframe that GYRE has to use
    npg_min, npg_max: int
        range in npg values to calculate modes
    gyre_base_file: string
        path to the template of basis gyre inlist to read and modify
    azimuthal_order, degree: int
        azimuthal order and degree of the modes
    """

    if gyre_summary_file == '':
        path, filename = gyre_in_file.rsplit('/',1)
        gyre_summary_file = f'{filename[:-3]}.HDF'

    with open(gyre_base_file, 'r') as f:
            lines = f.readlines()
    replacements = {
                    'FILENAME' : f'\'{mesa_pulsation_file}\'',
                    'OUTPUT'   : f'\'{gyre_summary_file}\'',
                    'N_PG_MIN' : f'{npg_min}',
                    'N_PG_MAX' : f'{npg_max}',
                    'FREQ_MIN' : f'{round(max(freq_min,0.01), 6)}',
                    'FREQ_MAX' : f'{round(freq_max, 6)}',
                    'GRIDFRAME': f'\'{rotation_frame}\'',
                    'OMEGAROT' : f'{omega_rot}',
                    'OMEGAUNIT': f'\'{unit_rot}\'',
                    'ORDER'    : f'{azimuthal_order}',
                    'DEGREE'   : f'{degree}'
                    }

    new_lines = []
    for line in lines:
            new_line = line
            for key in replacements:
                    if (replacements[key] != ''):
                            new_line = new_line.replace(key, replacements[key])
            new_lines.append(new_line)

    with open(gyre_in_file, 'w') as f:
            f.writelines(new_lines)
    return

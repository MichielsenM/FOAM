# Functions for building MESA and GYRE grids on the VSC framework
# from PyPulse import grid_building_vsc as gbv
import numpy as np
import glob, os, sys, csv
from . import my_python_functions as mypy
from . import functions_for_gyre as ffg
from shutil import copyfile

import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
logger = logging.getLogger('logger')
logger.setLevel(logging.DEBUG)

################################################################################
def make_mesa_setup(setup_directory=f'{os.getcwd()}/MESA_setup', work_dir=f'{os.getcwd()}/MESA_work_dir',
                    Z_ini_list=[0.014], M_ini_list=[1], log_Dmix_list=[1], aov_list=[0], fov_list=[0],
                    output_dir= os.path.expandvars(f'{os.getcwd()}/MESA_out')):
    """
    Construct a setup for a MESA grid with job lists to run on e.g. SLURM, and bash scripts to run each job list.
    ------- Parameters -------
    setup_directory, output_dir, work_dir: string
        paths to the directory with the MESA job submission files, the MESA output folder, and to the MESA work directory.
    Z_ini_list, M_ini_list, log_Dmix_list, fov_list, aov_list: list of floats
        Lists of the parameter values to be computed in the grid.
    """
    if 'site_scratch' in setup_directory:
        setup_directory = setup_directory[setup_directory.rfind("site_scratch"):].replace('site_scratch', '/scratch')
    if 'site_scratch' in work_dir:
        work_dir = work_dir[work_dir.rfind("site_scratch"):].replace('site_scratch', '/scratch')
    if 'site_scratch' in setup_directory:
        output_dir = output_dir[output_dir.rfind("site_scratch"):].replace('site_scratch', '/scratch')

    if not os.path.exists(work_dir):
        logger.error(f'Specified MESA work directory does not exist: {work_dir}')
        sys.exit()
    if not os.path.exists(f'{setup_directory}'):
        os.makedirs(f'{setup_directory}')

    with open(f'{setup_directory}/MESA_parameters.csv', 'w') as tsvfile:
        writer = csv.writer(tsvfile)
        header = ['Zini', 'Mini', 'logD', 'aov', 'fov', 'output_dir', 'MESA_work_dir']
        writer.writerow(header)

        for Z_ini in Z_ini_list:
            Z_ini = '{:5.3f}'.format(Z_ini)
            for M_ini in M_ini_list:
                M_ini = '{:4.2f}'.format(M_ini)
                for log_Dmix in log_Dmix_list:
                    log_Dmix ='{:4.2f}'.format(log_Dmix)
                    for a_ov in aov_list:
                        a_ov = '{:5.3f}'.format(a_ov)
                        for f_ov in fov_list:
                            f_ov = '{:5.3f}'.format(f_ov)

                            line_to_write = [Z_ini, M_ini, log_Dmix, a_ov, f_ov, output_dir, work_dir]
                            writer.writerow(line_to_write)

    copyfile(os.path.expandvars('$CONDA_PREFIX/lib/python3.7/site-packages/PyPulse/templates/run_MESA.sh'), f'{setup_directory}/run_MESA.sh')
    copyfile(os.path.expandvars('$CONDA_PREFIX/lib/python3.7/site-packages/PyPulse/templates/VSC_submit_MESA.pbs'), f'{setup_directory}/submit_MESA.pbs')
    return

################################################################################
def make_gyre_setup(setup_directory=f'{os.getcwd()}/GYRE_setup', npg_min=-50, npg_max=-1, azimuthal_order=1, degree=1, omega_rot=[0.0], unit_rot = 'CYC_PER_DAY', rotation_frame='INERTIAL',
                    output_dir=os.path.expandvars('$VSC_SCRATCH/GYRE_out'), mesa_dir=os.path.expandvars('$VSC_SCRATCH/MESA_out')):
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
    if not os.path.exists(f'{setup_directory}'):
        os.makedirs(f'{setup_directory}/inlists')

    gyre_files = glob.glob(mesa_dir + '/*/gyre/*.GYRE' )

    with open(f'{setup_directory}/GYRE_parameters.csv', 'w') as tsvfile:
        writer = csv.writer(tsvfile)
        header = ['Zini', 'Mini', 'logD', 'aov', 'fov', 'Xc', 'GYRE_inlist', 'output_dir']
        writer.writerow(header)

        for file_path in gyre_files:
            path, filename = file_path.rsplit('/',1)
            param_dict = mypy.get_param_from_filename(file_path, ['M', 'Z', 'logD', 'aov', 'fov', 'Xc'])
            # output_dir_Z = f'{output_dir}/Zini{param_dict["Z"]}/{filename[:filename.rfind(".")]}'
            for rotation in omega_rot:
                output_dir_Z = f'{output_dir}/rot{rotation}/Zini{param_dict["Z"]}'

                f_min, f_max = ffg.calc_scanning_range(file_path, npg_min=npg_min, npg_max=npg_max, l=degree, m=azimuthal_order, omega_rot=rotation, unit_rot=unit_rot, rotation_frame=rotation_frame)
                inlist_to_write = f'{setup_directory}/inlists/rot{rotation}_{filename[:-5]}.in'
                write_gyre_inlist(inlist_to_write, file_path, npg_min=npg_min,npg_max=npg_max, freq_min=f_min, freq_max=f_max, omega_rot=rotation, unit_rot=unit_rot, rotation_frame=rotation_frame)

                line_to_write = [param_dict["Z"], param_dict["M"], param_dict["logD"], param_dict["aov"], param_dict["fov"], param_dict["Xc"], inlist_to_write, output_dir_Z]
                writer.writerow(line_to_write)

    copyfile(os.path.expandvars('$CONDA_PREFIX/lib/python3.7/site-packages/PyPulse/templates/run_GYRE.sh'), f'{setup_directory}/run_GYRE.sh')
    copyfile(os.path.expandvars('$CONDA_PREFIX/lib/python3.7/site-packages/PyPulse/templates/VSC_submit_GYRE.pbs'), f'{setup_directory}/submit_GYRE.pbs')
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
                    'FILENAME' : '{}'.format("'"+mesa_pulsation_file+"'"),
                    'OUTPUT'   : '{}'.format("'"+gyre_summary_file+"'"),
                    'N_PG_MIN' : '{:1.0f}'.format(npg_min),
                    'N_PG_MAX' : '{:1.0f}'.format(npg_max),
                    'FREQ_MIN' : '{:8.6f}'.format(max(freq_min,0.01)),
                    'FREQ_MAX' : '{:8.6f}'.format(freq_max),
                    'GRIDFRAME': '{}'.format("'"+rotation_frame+"'"),
                    'OMEGAROT' : '{}'.format(float(omega_rot)),
                    'OMEGAUNIT': '{}'.format("'"+unit_rot+"'"),
                    'ORDER'    : '{:1.0f}'.format(azimuthal_order),
                    'DEGREE'   : '{:1.0f}'.format(degree)
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

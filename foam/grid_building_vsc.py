"""Functions for building MESA and GYRE grids on the VSC (Vlaams Supercomputer Centrum) framework."""
# from foam import grid_building_vsc as gbv
import numpy as np
import glob, os, sys, csv
import logging, pkgutil, multiprocessing
from shutil import copyfile
from pathlib import Path
from functools import partial
from foam import support_functions as sf
from foam import functions_for_gyre as ffg

logger = logging.getLogger('logger.gbv')
################################################################################
def make_mesa_setup(setup_directory=f'{os.getcwd()}/MESA_setup', work_dir=f'{os.getcwd()}/MESA_work_dir',
                    Z_ini_list=[0.014], M_ini_list=[1], log_Dmix_list=[1], aov_list=[0], fov_list=[0],
                    output_dir= os.path.expandvars(f'{os.getcwd()}/MESA_out')):
    """
    Construct a setup and job list to run a MESA grid on the VSC.
    ------- Parameters -------
    setup_directory, output_dir, work_dir: string
        paths to the directory with the MESA job submission files, the MESA output folder, and to the MESA work directory.
    Z_ini_list, M_ini_list, log_Dmix_list, fov_list, aov_list: list of floats
        Lists of the parameter values to be computed in the grid.
    """
    for directory_name in [setup_directory, work_dir, output_dir]:
        if 'site_scratch' in directory_name:
            directory_name = directory_name[directory_name.rfind("site_scratch"):].replace('site_scratch', '/scratch')

    if not os.path.exists(work_dir):
        logger.error(f'Specified MESA work directory does not exist: {work_dir}')
        sys.exit()
    Path(f'{setup_directory}').mkdir(parents=True, exist_ok=True)

    with open(f'{setup_directory}/MESA_parameters.csv', 'w') as tsvfile:
        writer = csv.writer(tsvfile)
        header = ['Zini', 'Mini', 'logD', 'aov', 'fov', 'output_dir', 'MESA_work_dir']
        writer.writerow(header)

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

                            line_to_write = [Z_ini, M_ini, log_Dmix, a_ov, f_ov, output_dir, work_dir]
                            writer.writerow(line_to_write)

    copyfile(os.path.expandvars(f'{Path(__file__).parent}/templates/run_MESA.sh'), f'{setup_directory}/run_MESA.sh')
    copyfile(os.path.expandvars(f'{Path(__file__).parent}/templates/VSC_submit_MESA.pbs'), f'{setup_directory}/submit_MESA.pbs')
    return

################################################################################
def make_gyre_setup(setup_directory=f'{os.getcwd()}/GYRE_setup', npg_min=-50, npg_max=-1, azimuthal_order=1, degree=1,
                    omega_rot=[0.0], unit_rot = 'CYC_PER_DAY', rotation_frame='INERTIAL',
                    output_dir=os.path.expandvars(f'{os.getcwd()}/GYRE_out'), mesa_dir=os.path.expandvars(f'{os.getcwd()}/MESA_out'),
                    gyre_base_inlist = None):
    """
    Construct a setup and job list to run a GYRE grid on the VSC.
    GYRE inlists and jobs will be created for each pulsation file found in the MESA directory.
    Scanning ranges in GYRE inlists will be set based on desired n_pg range and Asymptotic_dP in the MESA hist file.
    ------- Parameters -------
    setup_directory, output_dir, mesa_dir: string
        paths to the directory where the bash setup is being made,
        to the directory where the GYRE output will be stored, and to the MESA output directory.
    npg_min, npg_max, degree, azimuthal_order : int
        Quantum numbers of the modes to be calculated. Minimum and maximum radial order, degree (l) and azimuthal order (m)
    omega_rot: list of float
        rotation frequency of the model
    unit_rot: string
        unit of the rotation frequency, can be CYC_PER_DAY or CRITICAL (roche critical)
    rotation_frame: string
        rotational frame of reference for the pulsation freqencies
    gyre_base_inlist: string
        path to the template of basis gyre inlist to read and later modify, uses the template in this package by default
    """
    for directory_name in [setup_directory, mesa_dir, output_dir]:
        if 'site_scratch' in directory_name:
            directory_name = directory_name[directory_name.rfind("site_scratch"):].replace('site_scratch', '/scratch')

    if not os.path.exists(mesa_dir):
        logger.error(f'Specified MESA output directory does not exist: {mesa_dir}')
        sys.exit()
    Path(f'{setup_directory}/inlists').mkdir(parents=True, exist_ok=True)

    gyre_files = glob.glob(mesa_dir + '/*/gyre/*.GYRE' )

    if gyre_base_inlist is None:
        data = pkgutil.get_data(__name__, "templates/gyre_template.in")
        gyre_base_inlist_lines = data.decode(sys.stdout.encoding)
    else:
        with open(gyre_base_inlist, 'r') as f:
            gyre_base_inlist_lines = f.readlines()

    with open(f'{setup_directory}/GYRE_parameters.csv', 'w') as tsvfile:
        writer = csv.writer(tsvfile)
        header = ['Zini', 'Mini', 'logD', 'aov', 'fov', 'Xc', 'GYRE_inlist', 'output_dir']
        writer.writerow(header)

        p = multiprocessing.Pool()	# Multiprocessing pool, uses #processes = #CPUs
        for rotation in omega_rot:
            # set part of the input parameters already in 'func', and only iterate over 'gyre_files' in the 'p.imap'
            func = partial(gyre_process, output_dir=output_dir, setup_directory=setup_directory, npg_min=npg_min, npg_max=npg_max, degree=degree,
                           azimuthal_order=azimuthal_order, rotation=rotation, unit_rot=unit_rot, rotation_frame=rotation_frame, gyre_base_inlist_lines=gyre_base_inlist_lines)
            for result in p.imap(func, gyre_files):
                writer.writerow(result)

    copyfile(os.path.expandvars(f'{Path(__file__).parent}/templates/run_GYRE.sh'), f'{setup_directory}/run_GYRE.sh')
    copyfile(os.path.expandvars(f'{Path(__file__).parent}/templates/VSC_submit_GYRE.pbs'), f'{setup_directory}/submit_GYRE.pbs')
    return

################################################################################
def gyre_process(file_path, output_dir='', setup_directory='', npg_min=-50,npg_max=-1, degree=1, azimuthal_order=1,
                 rotation=0, unit_rot='CYC_PER_DAY', rotation_frame='INERTIAL', gyre_base_inlist_lines=None):
    """
    Write GYRE inlist and make a corresponding line for the CSV file with all the run parameters for submitting the jobs to the VSC
    ------- Parameters -------
    file_path: String
        path to the .GYRE file to construct an inlist for
    setup_directory, output_dir: string
        paths to the directory where the bash setup is being made, and the directory where the GYRE output will be stored.
    npg_min, npg_max, degree, azimuthal_order : int
        Quantum numbers of the modes to be calculated. Minimum and maximum radial order, degree (l) and azimuthal order (m)
    rotation: float
        rotation frequency of the model
    unit_rot: string
        unit of the rotation frequency, can be CYC_PER_DAY or CRITICAL (roche critical)
    rotation_frame: string
        rotational frame of reference for the pulsation freqencies
    gyre_base_inlist_lines: string
        lines of the template gyre inlist to modify

    ------- Returns -------
    line_to_write: string
        line to write in the CSV file containing all run parameters for submitting to the VSC
    """
    param_dict = sf.get_param_from_filename(file_path, ['M', 'Z', 'logD', 'aov', 'fov', 'Xc'])
    output_dir_Z = f'{output_dir}/rot{rotation}/Zini{param_dict["Z"]}_Mini{param_dict["M"]}'

    f_min, f_max = ffg.calc_scanning_range(file_path, npg_min=npg_min, npg_max=npg_max, l=degree, m=azimuthal_order, omega_rot=rotation,
                                           unit_rot=unit_rot, rotation_frame=rotation_frame)
    inlist_to_write = f'{setup_directory}/inlists/rot{rotation}_{Path(file_path).stem}.in'
    write_gyre_inlist(inlist_to_write, file_path, npg_min=npg_min,npg_max=npg_max, freq_min=f_min, freq_max=f_max, omega_rot=rotation,
                      unit_rot=unit_rot, rotation_frame=rotation_frame, gyre_base_inlist_lines=gyre_base_inlist_lines)

    line_to_write = [param_dict["Z"], param_dict["M"], param_dict["logD"], param_dict["aov"], param_dict["fov"], param_dict["Xc"], inlist_to_write, output_dir_Z]

    return line_to_write

################################################################################
def write_gyre_inlist( gyre_in_file, mesa_pulsation_file, gyre_summary_file='',
                       freq_min=0.01, freq_max=10, rotation_frame='INERTIAL',
                       npg_min=-50,npg_max=-1, omega_rot=0.0, unit_rot = 'CYC_PER_DAY', azimuthal_order=1, degree=1, gyre_base_inlist_lines = None,
                       gyre_base_inlist = None
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
    gyre_base_inlist_lines: string
        lines of the template gyre inlist to modify
    gyre_base_inlist: string
        used if gyre_base_inlist_lines was None,
        path to the template of basis gyre inlist to read and modify, uses the template in this package by default
    azimuthal_order, degree: int
        azimuthal order and degree of the modes
    """

    if gyre_summary_file == '':
        gyre_summary_file = f'{Path(gyre_in_file).stem}.HDF'

    if gyre_base_inlist_lines is None:
        if gyre_base_inlist is None:
            data = pkgutil.get_data(__name__, "templates/gyre_template.in")
            gyre_base_inlist_lines = data.decode(sys.stdout.encoding)
        else:
            with open(gyre_base_inlist, 'r') as f:
                gyre_base_inlist_lines = f.readlines()

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

    for key in replacements:
        if (replacements[key] != ''):
            gyre_base_inlist_lines = gyre_base_inlist_lines.replace(key, replacements[key])

    with open(gyre_in_file, 'w') as f:
        f.writelines(gyre_base_inlist_lines)

    return

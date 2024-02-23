""" Create a summary file to store e.g. all history files of a MESA grid in a nested dictionary."""
import sys, os
import logging
import hdfdict, h5py
import numpy as np

from foam import functions_for_mesa as ffm
from foam import support_functions as sf

logger = logging.getLogger('logger.mg')
################################################################################
def _make_nested_dict(list_keys, value):
    """
    Recursively make a nested dictionary with the keys at different levels set by list_keys,
    and the most nested level containing the value.
    
    Parameters
    ----------
    list_keys: list of keys
        List of keys for the different levels of the nested dictionary
    value: any
        The value, list, dictionary... that is coupled to the key

    Returns
    ----------
    dict: dict
        A nested dictionary
    """

    if len(list_keys) == 1:
        return {list_keys[0]: value}
    return {list_keys[0] : _make_nested_dict( list_keys[1:], value) }
################################################################################

class GridSummary:
    """
    Class to create and contain a summary of a MESA grid in a nested dictionary.
    """
    def __init__(self, grid_parameters=None, Zsun=0.014):
        """
        Summary of the information of the grid of stellar models.
        
        Parameters
        ----------
        grid_path: string
            Path to the stellar model grid.
        Zsun: float
            Solar metallicity value. Default = 0.014.
        """
        self.grid_parameters = grid_parameters
        self.grid_data = None
        self.Zsun = Zsun

    ############################################################################
    def create_summary_file(self, dir_path, columns=['star_age','log_L','log_R','log_Teff','log_g'],
            magnitudes=False, output_name='stellar_grid.h5', file_ending='hist', files_directory_name='history'):
        '''
        Create a hdf5 summary file containing certain columns of all history files in a MESA grid.
        
        Parameters
        ----------
        dir_path: string
            The path to the grid directory.
        columns: list of strings
            Names of the columns to be kept in the summary file.
        magnitudes: boolean
            If True, keep the absolute magnitude columns in the summary file.
        output_name: String
            Name for the summary file.
        '''
        if os.path.isfile(output_name):
            logger.warning(f'Output file exists already: {output_name}')
            sys.exit()

        summary_data ={}
        for subdir, dirs, files in os.walk(dir_path):
            # Skip the directory if there are no directories in it, and it's not the directory type we're looking for.
            if (len(dirs) == 0) and not subdir.endswith(files_directory_name):
                continue
            for file in files:
                if file.endswith(file_ending) and subdir.endswith(files_directory_name):
                    # Read in the files and get the main parameters from the filename
                    _, data = ffm.read_mesa_file(os.path.join(subdir, file))
                    filename_params = sf.get_param_from_filename(file, self.grid_parameters)

                    # Keep only the specified columns and the absolute magnitudes (if magnitudes==True)
                    for key in list(data.keys()):
                        if magnitudes == True and key.startswith('abs_mag_'):
                            continue
                        elif key not in columns:
                            del data[key]

                    # Create a nested dictionary with the order of levels the same as the grid_parameters
                    # and the column data dictionary at the bottom level
                    keys = [ filename_params[x] for x in self.grid_parameters]
                    for i in range(len(self.grid_parameters)):
                        param_list = keys
                        param = param_list.pop(0)

                        if i == 0:
                            if param not in summary_data:
                                summary_data[param] = _make_nested_dict(param_list, data)
                                break
                            else:
                                nested_summary = summary_data[param]

                        elif i == len(self.grid_parameters)-1:
                            nested_summary[param] = data

                        else:
                            if param not in nested_summary:
                                nested_summary[param] = _make_nested_dict(param_list, data)
                                break
                            else:
                                nested_summary = nested_summary[param]

        self.grid_data = summary_data
        self._set_param_ranges()
        grid_info = {'grid_data': summary_data , 'grid_parameters': self.grid_parameters}

        with h5py.File(  output_name, 'w' ) as hdf5_file:
            hdfdict.dump(grid_info, output_name)


    ############################################################################
    def read_summary_file(self, grid_path):
        '''
        Read in a stellar model grid in hdf5 format created with create_summary_file.

        Parameters
        ----------
        grid_path: string
            Path to the stellar model grid.
        '''
        # Load in the grid.
        grid = hdfdict.load(grid_path, lazy=False)
        self.grid_data = grid['grid_data']
        self.grid_parameters = [x.decode(sys.stdout.encoding) for x in grid['grid_parameters']]
        self._set_param_ranges()

    ############################################################################
    def _set_param_ranges(self):
        """
        Set the attributes holding the ranges of all grid parameters
        """
        dictionary = self.grid_data
        for parameter in self.grid_parameters:
            dict_keys = sorted(list(dictionary.keys()))
            values = np.array([x for x in dict_keys])
            dictionary = dictionary[dict_keys[-1]] # Highest mass values have largest range in logDext, hence use highest key of each parameter
            setattr(self, f'{parameter}_array', values)

"""
  This module provides different tools to read in the data form MESA consistently.
  It supports reading a single history log file, a single profile, or a large number of
  files from a repository.

  log history:
  - 19 August 2014: published on bitbucket
  - 25 July 2013: created.

"""
import sys,os
import logging
import numpy as np

#=================================================================================================================
def read_multiple_mesa_files(files, just_header=None, is_hist = False, is_prof = False):
  """
  'Reads "multiple" MESA hist/prof files. Returns a list of dictionaries.'
  'Each item in the list is a dictionary and each dictionary is associated to'
  'each single hist/prof file.'
  'Keys are the header and column titles.'
  'Values for header keys are floats, and values for columns are numpy arrays.'
  """
  list_flags = [is_hist, is_prof]
  if all(list_flags) or not any(list_flags):
    message = 'Error: read: read_multiple_mesa_files: set ONLY "is_hist" OR "is_prof" flag to True'
    raise SystemExit(message)
  if type(files) is not type([]):
    message = 'Error: read: read_multiple_mesa_files: You must provide a "list" of input filename(s)'
    raise SystemExit(message)

  n_files = len(files)
  if n_files == 0:
    print('File list is empty. Exit')
    raise SystemExit(0)

  logging.debug(' - mesa_gyre.read.read_multiple_mesa_files:')

  # read the first file, and extract the datatype, and pass the dtype to read_mesa_ascii
  header, data = read_mesa_ascii(files[0])
  dtype = data.dtype

  data_lis = []
  n_skip = 0
  for i_file, mesa_file in enumerate(files):
    dic = {}
    if not os.path.isfile(mesa_file):     # Insert an empty dictionary, to keep track of the file
      print(('   %s does not exist.' % (mesa_file, )))
      n_skip += 1
      dic['filename'] = ''
      dic['header'] = 0.0
      if is_hist: dic['hist'] = np.rec.fromarrays([0], names=['empty'])
      if is_prof: dic['prof'] = np.rec.fromarrays([0], names=['empty'])
      continue
    header, data = read_mesa_ascii(mesa_file, dtype=dtype)
    dic['filename'] = mesa_file
    dic['header'] = header
    if is_hist: dic['hist'] = data
    if is_prof: dic['prof'] = data
    data_lis.append(dic)
    if n_files<10:
      logging.debug(('   Read: %s, N_columns: %s, N_lines: %s' % (mesa_file, len(data.dtype.names), len(data))))

  if n_skip>0: print(('   %s files were skipped reading. \n' % (n_skip, )))
  else: logging.debug(('   All %s files successfully read.\n' % (n_files, )))

  return data_lis

#=================================================================================================================
def get_dtype(list_names):
  """
  Return hard-coded data types for data read from MESA history files or profiles.
  We treat integer columns as integer (np.int32), and all other float columns as single precision (np.float32),
  except if they belong to the list of double precision columns (np.float64).
  @param list_names: list of column names from the data block of the MESA history or profile file
  @type list_names: list of strings
  @return: numpy-compatible list of tuples for data type, each tuple is associated with a column, and has two elements.
       The first element is the column name, and the second is the numpy data type
  @rtype: list of tuples
  """
  n_col = len(list_names)
  if n_col == 0:
    message = 'get_dtype: empty input list of names'
    logging.error(message)
    raise SystemExit(message)

  set_int8_columns  = get_mesa_int8_col_names()    # for byte integer columns -128 < int < 128
  set_int16_columns = get_mesa_int16_col_names()   # for integer columns -32768 < int < 32767
  set_int32_columns = get_mesa_int32_col_names()   # for long integers -2147483648 < int < 2147483647
  set_dp_columns    = get_mesa_dp_col_names()

  list_tup_out = []
  for i_tup, name in enumerate(list_names):
    if name in set_int8_columns:
      list_tup_out.append( (name, np.int8) )
    elif name in set_int16_columns:
      list_tup_out.append( (name, np.int16) )
    elif name in set_int32_columns:
      list_tup_out.append( (name, np.int32) )
    elif name in set_dp_columns:
      list_tup_out.append( (name, np.float64) )
    else:
      list_tup_out.append( (name, np.float64) )

  return np.dtype(list_tup_out)

#=================================================================================================================
def get_mesa_int8_col_names():
  """
  Return the names of those columns that contain integer values in MESA history or profile files. They are listed
  in <mesa>/star/defaults/history_columns.list and <mesa>/star/defaults/profile_columns.list
  @return: set of strings with those columns that contain integer values in MESA
  @rtype: set
  """
  output = set([ 'mix_type_1', 'mix_type_2', 'mix_type_3', 'mix_type_4', 'mix_type_5', 'mix_type_6',             # hist
               'burn_type_1', 'burn_type_2', 'burn_type_3', 'burn_type_4', 'burn_type_5', 'burn_type_6',         # hist
               'mixing_type', 'mlt_mixing_type', 'sch_stable', 'ledoux_stable', 'stability_type'                 # prof
                ] )

  return output

#=================================================================================================================
def get_mesa_int16_col_names():
  """
  Return the names of those columns that contain integer values in MESA history or profile files. They are listed
  in <mesa>/star/defaults/history_columns.list and <mesa>/star/defaults/profile_columns.list
  @return: set of strings with those columns that contain integer values in MESA
  @rtype: set
  """
  output = set([ 'num_zones', 'cz_zone', 'cz_top_zone', 'num_backups', 'num_retries', 'zone'] )

  return output

#=================================================================================================================
def get_mesa_int32_col_names():
  """
  Return the names of those columns that contain integer values in MESA history or profile files. They are listed
  in <mesa>/star/defaults/history_columns.list and <mesa>/star/defaults/profile_columns.list
  @return: set of strings with those columns that contain integer values in MESA
  @rtype: set
  """
  output = set([ 'model_number', 'version_number', 'nse_fraction' ] )

  return output

#=================================================================================================================
def get_mesa_dp_col_names():
  """
  Return the names of those columns that contain double precision floats in MESA history or profile files. They are listed
  in <mesa>/star/defaults/history_columns.list and <mesa>/star/defaults/profile_columns.list
  These are those specific columns in the profiles that are used for the purpose of generating GYRE or GraCo input files.
  Since these codes require double precision input, these columns have to be preserved in double precision format.
  @return: set of strings with those columns that contain float values in MESA
  @rtype: set
  """
  output = set(['radius', 'mass', 'dq', 'eps_grav', 'cp', 'cv', 'Cv', 'brunt_N2_composition_term', 'lamb_S2',
                'q_div_xq', 'chiT_div_chiRho', 'kappa', 'kappa_rho', 'kappa_T', 'epsilon', 'epsilon_rho', 'epsilon_T',
                'omega', 'pressure', 'logP', 'pgas', 'prad', 'density', 'logRho', 'temperature', 'logT', 'luminosity',
                'luminosity_rad', 'luminosity_conv', 'gamma1', 'gradT', 'grada', 'brunt_N2', 'lnq', 'Cp', 'ln_free_e',
                'A2', 'Pg_div_P', 'gradr', 'h1', 'h2', 'he3', 'he4', 'li7', 'be7', 'c12', 'c13', 'n14', 'n15', 'o16',
                'o17', 'be9', 'si28' ])

  return output

#=================================================================================================================
def read_mesa_ascii(filename, dtype=None):
  """
  Read an history or profile ascii output from MESA.
  @param filename: full path to the input ascii file
  @type filename: string
  @param dtype: numpy-compatible dtype object. if it is not provided, it will be retrieved from read.get_dtype()
  @type dtype: list of tuples
  @return dictionary of the header of the file, and the record array for the data block. It can be called in the following way
     >>> header, data = read.read_mesa_ascii('filename')
  @rtype: dictionary and numpy record array
  """
  if not os.path.isfile(filename):
    message = 'read_mesa_ascii: {0} does not exist'.format(filename)
    logging.error(message)
    raise SystemExit(message)

  with open(filename, 'r') as r: lines = r.readlines()
  logging.info('read_mesa_ascii: {0} is read into list of lines'.format(filename))

  skip          = lines.pop(0)
  header_names  = lines.pop(0).rstrip('\r\n').split()
  header_strs   = lines.pop(0).rstrip('\r\n').split()
  header_vals   = []
  for val in header_strs:
    if 'D' in val:
      val = val.replace('D', 'E')
    if '"' in val:
      val = 0
    header_vals.append(val)
  temp          = np.array([header_vals], float).T
  header        = np.core.records.fromarrays(temp, names=header_names)
  skip          = lines.pop(0)
  skip          = lines.pop(0)

  col_names     = lines.pop(0).rstrip('\r\n').split()
  n_cols        = len(col_names)
  if dtype is None:
    col_dtype     = get_dtype(col_names)
  else:
    col_dtype   = dtype
  if n_cols != len(col_dtype):
    message = 'read_mesa_ascii: Incompatible number of provided datatype objects'
    logging.error(message)
    raise SystemExit(message)

  data          = []
  for i_line, line in enumerate(lines):
    if not line.rstrip('\r\n').split(): continue  # skip empty lines
    line = line.replace('D', 'E')
    data.append(line.rstrip('\r\n').split())

  data = np.core.records.fromarrays(np.array(data, float).transpose(), dtype=col_dtype)

  return header, data
#=================================================================================================================

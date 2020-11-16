""" Top level script to run the pipeline sequentially,
comment certain imports if you don't want to repeat them on repeated runs."""
import importlib

importlib.import_module('PyPulse.pipeline.0_extract_puls&spectro')
importlib.import_module('PyPulse.pipeline.1_constuct_pattern')
importlib.import_module('PyPulse.pipeline.2_calculate_likelihood')
importlib.import_module('PyPulse.pipeline.3_spectroClip_AICc')
importlib.import_module('PyPulse.pipeline.4_bestModel_errors')
importlib.import_module('PyPulse.pipeline.5_correlationPlots')

importlib.import_module('PyPulse.pipeline.table_bestModels')

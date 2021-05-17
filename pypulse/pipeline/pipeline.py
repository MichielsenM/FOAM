""" Top level script to run the pipeline sequentially,
comment certain imports if you don't want to repeat them on repeated runs."""
import importlib

importlib.import_module('pypulse.pipeline.0_extract_puls&spectro')
importlib.import_module('pypulse.pipeline.1_constuct_pattern')
importlib.import_module('pypulse.pipeline.2_calculate_likelihood')
importlib.import_module('pypulse.pipeline.3_spectroClip_AICc')
importlib.import_module('pypulse.pipeline.4_bestModel_errors')
importlib.import_module('pypulse.pipeline.5_correlationPlots')

importlib.import_module('pypulse.pipeline.table_bestModels')

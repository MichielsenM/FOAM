""" Top level script to run the pipeline sequentially, copy this script to the folder where you want to run the analysis.
Comment specific imports if you don't want to repeat them on repeated runs."""
import importlib, os
from pathlib import Path
from foam.pipeline import pipelineConfig

pipelineConfig.config = pipelineConfig.pipelineConfig()

# Run the pipeline
importlib.import_module('foam.pipeline.pipe0_extract_puls_and_spectro')
importlib.import_module('foam.pipeline.pipe1_construct_pattern')

# Change the current working directory for nested grids
if pipelineConfig.config.fixed_parameters is not None:
    Path(pipelineConfig.config.nested_grid_dir).mkdir(parents=True, exist_ok=True)
    os.chdir(pipelineConfig.config.nested_grid_dir)

importlib.import_module('foam.pipeline.pipe2_calculate_likelihood')
importlib.import_module('foam.pipeline.pipe3_spectroClip_AICc')
importlib.import_module('foam.pipeline.pipe4_bestModel_errors')
importlib.import_module('foam.pipeline.pipe5_correlationPlots')
importlib.import_module('foam.pipeline.pipe6_table_bestModels')

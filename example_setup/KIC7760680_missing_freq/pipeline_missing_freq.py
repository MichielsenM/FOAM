""" Top level script to run the pipeline sequentially, copy this script to the folder where you want to run the analysis.
Comment specific imports if you don't want to repeat them on repeated runs."""
if __name__ == '__main__':
    import importlib, os
    from pathlib import Path
    from foam.pipeline import pipeline_config

    pipeline_config.config = pipeline_config.PipelineConfig(star = 'KIC7760680',
                                    observations='data_KIC7760680_missing_freq.tsv',
                                    pattern_starting_pulsation = {'period': [0.90147, 0.981691, 1.094833, 1.46046], 'frequency': [None]},
                                    grid_parent_directory = None, #change to '/YOUR_PATH/MESA_grid' if you need to perform step 0
                                    grids = ['DO'], #Diffusive Overshooting grid
                                    rotation_gyre=0.4805,
                                    N_periods = 32,
                                    N_pattern_parts = 4,
                                    pattern_methods = ['provided-pulsation', 'highest-frequency'], # exclude chisq_longest_sequence for faster runtime
                                    observable_seismic = ['P'],
                                    nr_cpu = 4)

    # Run the pipeline
    pipeline_config.config.logger.info('step 0: Extracting grid')
    importlib.import_module('foam.pipeline.pipe0_extract_grid')
    pipeline_config.config.logger.info('step 0: Done\n')

    pipeline_config.config.logger.info('1: Constructing theoretical patterns')    
    importlib.import_module('foam.pipeline.pipe1_construct_pattern')
    pipeline_config.config.logger.info('step 1: Done\n')
    
    # Change the current working directory for nested grids
    if pipeline_config.config.fixed_parameters is not None:
        Path(pipeline_config.config.nested_grid_dir).mkdir(parents=True, exist_ok=True)
        os.chdir(pipeline_config.config.nested_grid_dir)
    
    pipeline_config.config.logger.info('2: Calculating liklihoods')    
    importlib.import_module('foam.pipeline.pipe2_calculate_likelihood')
    pipeline_config.config.logger.info('step 2: Done\n')

    pipeline_config.config.logger.info('3: Adding constraints')
    importlib.import_module('foam.pipeline.pipe3_add_constraints')
    pipeline_config.config.logger.info('step 3: Done\n')

    pipeline_config.config.logger.info('4: Calculating AIC')
    importlib.import_module('foam.pipeline.pipe4_AICc')
    pipeline_config.config.logger.info('step 4: Done\n')

    pipeline_config.config.logger.info('5: Calculating best model errors')
    importlib.import_module('foam.pipeline.pipe5_best_model_errors')
    pipeline_config.config.logger.info('step 5: Done\n')

    pipeline_config.config.logger.info('6: Backing plots into a corner')
    importlib.import_module('foam.pipeline.pipe6_corner_plots')
    pipeline_config.config.logger.info('step 6: Done\n')

    pipeline_config.config.logger.info('7: Making table with best models')
    importlib.import_module('foam.pipeline.pipe7_table_best_models')
    pipeline_config.config.logger.info('step 7: Done\n')
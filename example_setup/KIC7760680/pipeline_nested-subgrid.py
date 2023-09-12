""" Top level script to run the pipeline sequentially, copy this script to the folder where you want to run the analysis.
Comment specific imports if you don't want to repeat them on repeated runs."""
if __name__ == '__main__':
    import importlib, os
    from pathlib import Path
    from foam.pipeline import pipeline_config

    pipeline_config.config = pipeline_config.PipelineConfig(star = 'KIC7760680',
                                    observations='data_KIC7760680.tsv',
                                    highest_amplitude_pulsation = {'period': [1.158919], 'frequency': [None]},
                                    grid_parent_directory = '/YOUR_PATH/MESA_grid',
                                    grids = ['DO'], #Diffusive Overshooting grid
                                    rotation_gyre=0.4805,
                                    N_periods = 36,
                                    observable_seismic = ['dP'],
                                    free_parameters = ['Z', 'M', 'logD', 'Xc'],
                                    fixed_parameters = {'aov':0, 'fov':0},
                                    nr_cpu = 4)


    # Run the pipeline
    importlib.import_module('foam.pipeline.pipe0_extract_grid')
    print('step 0 done')
    importlib.import_module('foam.pipeline.pipe1_construct_pattern')
    print('step 1 done')
    
    # Change the current working directory for nested grids
    if pipeline_config.config.fixed_parameters is not None:
        Path(pipeline_config.config.nested_grid_dir).mkdir(parents=True, exist_ok=True)
        os.chdir(pipeline_config.config.nested_grid_dir)

    importlib.import_module('foam.pipeline.pipe2_calculate_likelihood')
    print('step 2 done')
    importlib.import_module('foam.pipeline.pipe3_add_constraints')
    print('step 3 done')
    importlib.import_module('foam.pipeline.pipe4_AICc')
    print('step 4 done')
    importlib.import_module('foam.pipeline.pipe5_best_model_errors')
    print('step 5 done')
    importlib.import_module('foam.pipeline.pipe6_corner_plots')
    print('step 6 done')
    importlib.import_module('foam.pipeline.pipe7_table_best_models')
    print('step 7 done')

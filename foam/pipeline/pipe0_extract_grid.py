""" From the grid location, extract all relevant info from MESA profiles and all frequencies from GYRE models."""

from pathlib import Path

from foam import functions_for_gyre as ffg
from foam import functions_for_mesa as ffm
from foam.pipeline.pipeline_config import config

################################################################################

for grid in config.grids:
    output_file = f"../grid_summary/surfaceGrid_{grid}.hdf"
    if not Path(output_file).is_file():
        ffm.extract_surface_grid(
            f"{config.grid_parent_directory}/{grid}/MESA_out/*/profiles/*{config.subgrid}*prof",
            output_file=output_file,
            parameters=config.grid_parameters,
            nr_cpu=config.nr_cpu,
            additional_observables=config.observable_additional,
        )
    else:
        config.logger.warning(f"file already existed: {output_file}")

    output_file = f"../grid_summary/pulsationGrid_{grid}_rot{config.rotation_gyre}_k{config.kval}m{config.mval}.hdf"
    if not Path(output_file).is_file():
        ffg.extract_frequency_grid(
            f"{config.grid_parent_directory}/{grid}/GYRE_out/rot{config.rotation_gyre}_k{config.kval}m{config.mval}/*/*{config.subgrid}*.HDF",
            output_file=output_file,
            parameters=["rot"] + config.grid_parameters,
            nr_cpu=config.nr_cpu,
        )
    else:
        config.logger.warning(f"file already existed: {output_file}")

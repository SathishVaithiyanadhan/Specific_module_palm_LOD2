# PALM Specific emission input

# Specific Emission Mode Input Data in LOD 2

This repository provides a modular workflow to generate the Specific Emission Mode Input Data in LOD 2 for the [PALM modeling system](https://gitlab.palm-model.org/releases/palm_model_system/-/releases) using downscaled GRETA emission inventories. The scripts process spatiotemporal emission data into CF-compliant NetCDF files compatible with PALM's LOD2 input for traffic emissions will be called [model]_emis_traffic, while that for the all-purpose generic mode will be called [model]_emis_generic.

LOD 2: Gridded preprocessed hourly (other temporal intervals will be possible in later versions) emission information that is already temporally disaggregated must be supplied by the user. IMPORTANT: In this mode, the initial date of the simulation has to coincide with the first day for which emission values are available - source: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/chemdesc

---

# Attributes and Dimensions

To setup the attributes and the dimensions for the specific emission driver, following the latest [Specific Emission Mode Input Data in LOD 2](https://docs.palm-model.com/23.04/Guide/LES_Model/Modules/Chemistry/EMISSIONS_LOD2_spec/) documentation. This emission input follows the PIDS for the PALM version 25.04. 

## Key features

The driver integrates gridded emission data (e.g., from the GRETA inventory) with PALM's urban microclimate simulations. It implements:
- **AOI extraction** and grids verification from the input static data (the static driver for the simulation is read, necessary data is extracted.)
- **Multiple sector handling**  from the GRETA emission inventory.
- **Multiple species handling**  from the GRETA emission inventory.
- **Hourly emissions** as the input (LOD2).
- **Dynamic unit conversion** includes all volumetric emission sources are expressed in mol/(m3s) for gas-phase species and kg/(m3s) for particulate matter (PM).
- Automatic **volumetric emission sources** for each species in the config, defining all source cell locations (i,j,k) defined in the variables vsrc_i, vsrc_j, and vsrc_k. 
- Properly handling the NAN values.

---

## Input data

The following data is required to create the chemistry driver for the PALM simulation using this tool.

1. Downscaled GRETA Emission inventory
	* Check the repo downscale_emissions_local **(https://git.rz.uni-augsburg.de/vaithisa/downscale_emissions_local.git)** to create your own input data. 

2. Static Data 
	* The static data which you have created using the Geospatial data to describe the topography of the simulation domain. 
    * It is used here to extract the AOI and Grid details for the PALM simulation. 

---

## Usage

1. Configure Paths/Parameters
   - Edit Specific_module_config.py:

       * Set **geotiff_dir**, **static_pth** and **output_dir**

       * Select the preferred **active categories** and **species**

       * Keep the **selected_band_prefix = "SumAllSectors"** to process all the active categories.

       * Set the **layer_height** from your PALM configuration.

       * Set the **non_zero_threshold** value to find the emission source locations.

       * Set the **model_name** (name of you PALM job) and **emission_mode** (Specific mode like traffic, generic, etc.,) to name the output file accordingly.

2. Run Main Script

    * **python Specific_module_tool.py** 

3. Output

    * NetCDF file generated at: **output_dir**

---

## Authors and acknowledgment

Show your appreciation to those who have contributed to the project.
For details and comments, please contact:
1. Sathish Kumar Vaithiyanadhan (sathish.vaithiyanadhan@uni-a.de)
2. Christoph Knote (christoph.knote@med.uni-augsburg.de)

@ Chair of Model-based Environmental Exposure Science (MBEES), Faculty of Medicine, University of Augsburg, Germany.

---

## License

For open source projects, say how it is licensed.

---

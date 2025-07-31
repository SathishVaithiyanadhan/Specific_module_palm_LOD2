# PALM Specific emission input

# Specific Emission Mode Input Data in LOD 2

This repository provides a modular workflow to generate the Specific Emission Mode Input Data in LOD 2 for the [PALM modeling system](https://gitlab.palm-model.org/releases/palm_model_system/-/releases) using downscaled GRETA emission inventories. The scripts process spatiotemporal emission data into CF-compliant NetCDF files compatible with PALM's LOD2 input for traffic emissions will be called [model]_emis_traffic, while that for the all-purpose generic mode will be called [model]_emis_generic.
---

# Attributes and Dimensions

To setup the attributes and the dimenions for the specific emission driver, follow the latest Specific Emission Mode Input Data in LOD 2 documentation in [PALM Specific emission input](https://docs.palm-model.com/23.04/Guide/LES_Model/Modules/Chemistry/EMISSIONS_LOD2_spec/). This emission input follows the PIDS for the PALM version 25.04. 

## Key features

The driver integrates gridded emission data (e.g., from the GRETA inventory) with PALM's urban microclimate simulations. It implements:
- **AOI extarction** and grids verfication from the input static data (the static driver for the simualtion is read, necessary data is extracted.)
- **Multiple sector handling**  from the GRETA emission inventory.
- **Multiple species handling**  from the GRETA emission inventory.
- **Hourly emissions** as the input (LOD2)
- **Dynamic unit conversion** includes all volumetric emission sources are expressed in mol/(m3s) for gas-phase species and kg/(m3s) for particulate matter (PM).
- Automatic **volumetric emission sources** for each species in the config, defining all source cell locations (i,j,k) defined in the varibles vsrc_i, vsrc_j, and vsrc_k. 
- Properly handling the NAN values.

---

## Input data

The following data is required to create the chemistry driver for the PALM simulation using this tool.

1. Downscaled GRETA Emission inventory
	* Check the repo downscale_emissions_local **(https://git.rz.uni-augsburg.de/vaithisa/downscale_emissions_local.git)** to create your own input data. 

2. Static Data 
	* The static data whcih you have created using the Geospatial data to describe the topography of the simulation domain. 
    * It is used here to extract the AOI and Grid details for the PALM simulation. 

---

## Usage

1. Configure Paths/Parameters
   - Edit Specific_module_config.py:

       * Set emis_geotiff_pth and static_pth 

       * Select the preferred active categories and species

2. Run Main Script

    * **python Specific_module_main.py** 

3. Output

    * NetCDF file generated at: **{static_pth}/{static}_emis_traffic**

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

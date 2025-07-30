import os
import re
import numpy as np
import netCDF4 as nc
from osgeo import gdal
from Specific_module_config import *
from Specific_module_main import resample_geotiff, read_geotiff_band
from collections import defaultdict
import datetime

def parse_band_description(desc):
    pattern = r'^([A-Za-z_]+)_h(\d+)_(\d{8})$'
    match = re.match(pattern, desc)
    return match.groups() if match else None

def get_time_bands_info(geotiff_path):
    ds = gdal.Open(geotiff_path)
    if not ds: raise ValueError(f"Can't open {geotiff_path}")
    
    time_info = defaultdict(lambda: defaultdict(list))
    for band_num in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_num)
        if desc := band.GetDescription():
            if (parsed := parse_band_description(desc)) and parsed[0] in active_categories:
                time_info[parsed[2]][f"h{parsed[1]}"].append({
                    'band_num': band_num,
                    'category': parsed[0],
                    'hour_num': int(parsed[1])
                })
    ds = None
    return time_info

def create_time_dimensions(time_info):
    dates = sorted(time_info.keys())
    hours = sorted({h for d in dates for h in time_info[d]}, key=lambda x: int(x[1:]))
    return [{'date': d, 'hour': h, 'hour_num': int(h[1:])} 
            for d in dates for h in hours if h in time_info[d]]

def get_source_locations(static_params, emission_data, spec_name_str):
    """Identify source cell locations (i,j,k) where emissions or building heights are non-zero."""
    ny, nx = static_params['ny'], static_params['nx']
    source_locations = []
    
    building_height = static_params['building_height']
    for j in range(ny):
        for i in range(nx):
            has_emission = False
            for spec in spec_name_str:
                for ts_data in emission_data.get(spec, {}).values():
                    for band_data in ts_data.values():
                        for band in band_data:
                            arr = read_geotiff_band(
                                f"{emis_geotiff_pth}emission_{spec}_temporal.tif",
                                band['band_num'],
                                static_params
                            )
                            if not np.isnan(arr[j, i]) and arr[j, i] != -9999.9 and arr[j, i] > 0:
                                has_emission = True
                                break
                        if has_emission:
                            break
                    if has_emission:
                        break
            if building_height[j, i] > 0 or has_emission:
                source_locations.append((i, j, 0))
    
    return source_locations

def convert_emission_units(arr, spec, dz):
    """Convert emission units from kg/m²/hour to mol/(m³s) for gases or kg/(m³s) for PM."""
    if spec in molar_mass and molar_mass[spec] is not None:
        arr_mol = arr / (molar_mass[spec] / 1000)
        arr_mol_per_m3s = arr_mol / (dz * 3600)
        return arr_mol_per_m3s
    else:
        return arr / (dz * 3600)

def create_chemistry_driver(static_params):
    print("\nCreating PALM LOD2 Chemistry Driver for Traffic Emissions")
    
    # Define species name mapping to uppercase format
    species_mapping = {spec: spec.upper().replace('_', '') for spec in spec_name_str}
    
    # Convert species names to uppercase format
    uppercase_spec_names = [species_mapping[spec] for spec in spec_name_str]
    
    # Collect temporal emission data for specified species
    all_time_info = {}
    for spec in spec_name_str:
        gt_path = f"{emis_geotiff_pth}emission_{spec}_temporal.tif"
        if os.path.exists(gt_path):
            all_time_info[spec] = get_time_bands_info(gt_path)
    
    if not all_time_info: 
        raise RuntimeError("No valid emission files found for specified species")

    # Get source locations (i,j,k)
    source_locations = get_source_locations(static_params, all_time_info, spec_name_str)
    if not source_locations:
        raise RuntimeError("No valid emission source locations found")
    
    # Create time dimensions
    time_steps = create_time_dimensions(next(iter(all_time_info.values())))
    
    # Create NetCDF dataset
    with nc.Dataset(f"{static_pth}{static}_emis_traffic", 'w') as ds:
        # Global attributes
        ds.setncatts({
            'description': 'LOD2 traffic emissions for PALM model simulation from Downscaled GRETA Emissions',
            'author': 'Sathish Kumar Vaithiyanadhan (sathish.vaithiyanadhan@med.uni-augsburg.de)',
            'Institution': 'Chair of Model-based Environmental Exposure Science (MBEES), Faculty of Medicine, University of Augsburg, Germany.',
            'lod': 2,
            'origin_time': static_params['origin_time'],
            'origin_lat': static_params['origin_lat'],
            'origin_lon': static_params['origin_lon'],
            'origin_x': static_params['origin_x'],
            'origin_y': static_params['origin_y'],
            'resolution': static_params['dx'],
            'history': f'Created on {datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S +000")}',
            'source': 'multiple sources; netCDF3',
            'origin_z': static_params['z_origin'],
            'Conventions': 'CF-1.7'
        })

        # Dimensions
        ds.createDimension('ntime', None)  # Unlimited
        ds.createDimension('nspecies', len(spec_name_str))
        ds.createDimension('nvsrc', len(source_locations))
        ds.createDimension('field_length', 64)

        # Coordinate variables
        vsrc_i = ds.createVariable('vsrc_i', 'i4', ('nvsrc',))
        vsrc_i.setncatts({'description': 'i grid indices for volume source location'})
        vsrc_i[:] = [loc[0] for loc in source_locations]

        vsrc_j = ds.createVariable('vsrc_j', 'i4', ('nvsrc',))
        vsrc_j.setncatts({'description': 'j grid indices for volume source location'})
        vsrc_j[:] = [loc[1] for loc in source_locations]

        vsrc_k = ds.createVariable('vsrc_k', 'i4', ('nvsrc',))
        vsrc_k.setncatts({'description': 'k grid indices for volume source location'})
        vsrc_k[:] = [loc[2] for loc in source_locations]

        # Time variables
        time = ds.createVariable('timestamp', 'S1', ('ntime', 'field_length'))
        time.setncatts({'description': 'Time stamps'})

        # Species variable
        species = ds.createVariable('species', 'S1', ('nspecies', 'field_length'))
        species.setncatts({'description': 'Emission species'})
        species[:] = nc.stringtochar(np.array(uppercase_spec_names, dtype='S64'))

        # Emission variables
        emission_vars = {}
        gas_phase_species = ['N2O', 'NOX', 'NMVOC', 'SO2', 'CO', 'NH3', 'CO2', 'CH4', 'NO', 'NO2', 'O3']
        for spec in uppercase_spec_names:
            var_name = f'vsrc_{spec}'
            emission_vars[spec] = ds.createVariable(var_name, 'f4', ('ntime', 'nvsrc'),
                                                  fill_value=np.float32(0.0))
            emission_vars[spec].setncatts({
                'description': f'volume source values for {spec}',
                'units': 'mol/(m^3 s)' if spec in gas_phase_species else 'kg/(m^3 s)',
                'missing_value': np.float32(0.0)
            })

        # Process emissions data
        for ts_idx, ts in enumerate(time_steps):
            dt = datetime.datetime.strptime(ts['date'], "%Y%m%d") + datetime.timedelta(hours=ts['hour_num'])
            formatted_timestamp = dt.strftime("%Y-%m-%d %H:%M:%S +000")
            time[ts_idx] = nc.stringtochar(np.array(formatted_timestamp.ljust(64), dtype='S64'))
            
            for spec_idx, (spec, upper_spec) in enumerate(zip(spec_name_str, uppercase_spec_names)):
                emission_array = np.zeros(len(source_locations), dtype=np.float32)
                if spec in all_time_info:
                    total_emission = np.zeros((static_params['ny'], static_params['nx']), 
                                            dtype=np.float32)
                    for band in all_time_info[spec][ts['date']][ts['hour']]:
                        arr = read_geotiff_band(
                            f"{emis_geotiff_pth}emission_{spec}_temporal.tif",
                            band['band_num'],
                            static_params
                        )
                        valid_mask = ~np.isnan(arr) & (arr != -9999.9)
                        total_emission[valid_mask] = np.where(
                            total_emission[valid_mask] == 0,
                            arr[valid_mask],
                            total_emission[valid_mask] + arr[valid_mask]
                        )
                    
                    total_emission = convert_emission_units(total_emission, spec, static_params['dz'])
                    
                    for src_idx, (i, j, k) in enumerate(source_locations):
                        emission_array[src_idx] = total_emission[j, i]
                
                emission_vars[upper_spec][ts_idx, :] = emission_array

    print("\nSuccessfully created PALM LOD2 traffic emissions driver with:")
    print(f"- {len(spec_name_str)} species")
    print(f"- {len(time_steps)} time steps")
    print(f"- {len(source_locations)} source locations")


###
#parallel procesing
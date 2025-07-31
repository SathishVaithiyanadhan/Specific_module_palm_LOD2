import os
import re
import numpy as np
import netCDF4 as nc
from osgeo import gdal
from Specific_module_config import *
from Specific_module_main import resample_geotiff, read_geotiff_band
from collections import defaultdict
import datetime
from multiprocessing import Pool, cpu_count
from functools import partial
import concurrent.futures

def parse_band_description(desc):
    """Compile regex pattern once for better performance"""
    pattern = re.compile(r'^([A-Za-z_]+)_h(\d+)_(\d{8})$')
    match = pattern.match(desc)
    return match.groups() if match else None

def get_time_bands_info(geotiff_path):
    """Optimized version with better memory management"""
    ds = gdal.Open(geotiff_path)
    if not ds:
        raise ValueError(f"Can't open {geotiff_path}")
    
    # Use regular dict with setdefault instead of lambda defaultdict
    time_info = {}
    for band_num in range(1, ds.RasterCount + 1):
        band = ds.GetRasterBand(band_num)
        if desc := band.GetDescription():
            if (parsed := parse_band_description(desc)) and parsed[0] in active_categories:
                date_dict = time_info.setdefault(parsed[2], {})
                hour_list = date_dict.setdefault(f"h{parsed[1]}", [])
                hour_list.append({
                    'band_num': band_num,
                    'category': parsed[0],
                    'hour_num': int(parsed[1])
                })
    del ds  # Explicit cleanup
    return time_info

def create_time_dimensions(time_info):
    """Optimized with list comprehensions"""
    dates = sorted(time_info.keys())
    hours = sorted({h for date in dates for h in time_info[date]}, key=lambda x: int(x[1:]))
    return [{'date': d, 'hour': h, 'hour_num': int(h[1:])} 
            for d in dates for h in hours if h in time_info[d]]

def process_cell(args):
    """Optimized cell processing for parallel execution"""
    i, j, static_params, emission_data, spec_name_str, emis_geotiff_pth = args
    building_height = static_params['building_height']
    
    # Early exit if building height > 0
    if building_height[j, i] > 0:
        return (i, j, 0)
    
    # Check emissions only if building height is 0
    for spec in spec_name_str:
        if spec not in emission_data:
            continue
        for date_data in emission_data[spec].values():
            for hour_data in date_data.values():
                for band in hour_data:
                    try:
                        arr = read_geotiff_band(
                            f"{emis_geotiff_pth}emission_{spec}_temporal.tif",
                            band['band_num'],
                            static_params
                        )
                        val = arr[j, i]
                        if not np.isnan(val) and val != -9999.9 and val > 0:
                            return (i, j, 0)
                    except Exception:
                        continue
    return None

def get_source_locations(static_params, emission_data, spec_name_str):
    """Parallel implementation with chunking for better performance"""
    ny, nx = static_params['ny'], static_params['nx']
    args = [(i, j, static_params, emission_data, spec_name_str, emis_geotiff_pth) 
            for j in range(ny) for i in range(nx)]
    
    # Use ProcessPoolExecutor for better resource management
    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        results = list(executor.map(process_cell, args, chunksize=max(1, len(args)//(cpu_count()*4))))
    
    return [loc for loc in results if loc is not None]

def convert_emission_units(arr, spec, dz):
    """Vectorized unit conversion"""
    if spec in molar_mass and molar_mass[spec] is not None:
        return arr / (molar_mass[spec] / 1000) / (dz * 3600)
    return arr / (dz * 3600)

def process_species_time_step(args):
    """Parallel processing of species for each time step"""
    ts, spec, upper_spec, all_time_info, static_params, source_locations = args
    emission_array = np.zeros(len(source_locations), dtype=np.float32)
    
    if spec in all_time_info and ts['date'] in all_time_info[spec] and ts['hour'] in all_time_info[spec][ts['date']]:
        total_emission = np.zeros((static_params['ny'], static_params['nx']), dtype=np.float32)
        
        for band in all_time_info[spec][ts['date']][ts['hour']]:
            try:
                arr = read_geotiff_band(
                    f"{emis_geotiff_pth}emission_{spec}_temporal.tif",
                    band['band_num'],
                    static_params
                )
                valid_mask = ~np.isnan(arr) & (arr != -9999.9)
                np.maximum(total_emission, arr, out=total_emission, where=valid_mask)
            except Exception:
                continue
        
        total_emission = convert_emission_units(total_emission, spec, static_params['dz'])
        
        # Vectorized assignment
        j_indices = [loc[1] for loc in source_locations]
        i_indices = [loc[0] for loc in source_locations]
        emission_array[:] = total_emission[j_indices, i_indices]
    
    return upper_spec, emission_array

def create_chemistry_driver(static_params):
    print("\nCreating PALM LOD2 Chemistry Driver for Traffic Emissions")
    
    # Precompute species mappings
    species_mapping = {spec: spec.upper().replace('_', '') for spec in spec_name_str}
    uppercase_spec_names = [species_mapping[spec] for spec in spec_name_str]
    gas_phase_species = {'N2O', 'NOX', 'NMVOC', 'SO2', 'CO', 'NH3', 'CO2', 'CH4', 'NO', 'NO2', 'O3'}
    
    # Parallel collection of temporal emission data
    all_time_info = {}
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_spec = {
            executor.submit(get_time_bands_info, f"{emis_geotiff_pth}emission_{spec}_temporal.tif"): spec
            for spec in spec_name_str if os.path.exists(f"{emis_geotiff_pth}emission_{spec}_temporal.tif")
        }
        for future in concurrent.futures.as_completed(future_to_spec):
            spec = future_to_spec[future]
            all_time_info[spec] = future.result()
    
    if not all_time_info:
        raise RuntimeError("No valid emission files found for specified species")

    # Get source locations in parallel
    source_locations = get_source_locations(static_params, all_time_info, spec_name_str)
    if not source_locations:
        raise RuntimeError("No valid emission source locations found")
    
    time_steps = create_time_dimensions(next(iter(all_time_info.values())))
    
    # Create NetCDF dataset
    with nc.Dataset(f"{static_pth}{static}_emis_traffic", 'w') as ds:
        # Set global attributes
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

        # Define dimensions
        ds.createDimension('ntime', None)
        ds.createDimension('nspecies', len(spec_name_str))
        ds.createDimension('nvsrc', len(source_locations))
        ds.createDimension('field_length', 64)

        # Create variables
        vsrc_i = ds.createVariable('vsrc_i', 'i4', ('nvsrc',))
        vsrc_i[:] = [loc[0] for loc in source_locations]
        
        vsrc_j = ds.createVariable('vsrc_j', 'i4', ('nvsrc',))
        vsrc_j[:] = [loc[1] for loc in source_locations]
        
        vsrc_k = ds.createVariable('vsrc_k', 'i4', ('nvsrc',))
        vsrc_k[:] = [loc[2] for loc in source_locations]
        
        time = ds.createVariable('timestamp', 'S1', ('ntime', 'field_length'))
        
        species = ds.createVariable('species', 'S1', ('nspecies', 'field_length'))
        species[:] = nc.stringtochar(np.array(uppercase_spec_names, dtype='S64'))

        # Create emission variables
        emission_vars = {}
        for spec in uppercase_spec_names:
            units = 'mol/(m^3 s)' if spec in gas_phase_species else 'kg/(m^3 s)'
            emission_vars[spec] = ds.createVariable(f'vsrc_{spec}', 'f4', ('ntime', 'nvsrc'),
                                                  fill_value=np.float32(0.0))
            emission_vars[spec].units = units

        # Process time steps in parallel
        for ts_idx, ts in enumerate(time_steps):
            dt = datetime.datetime.strptime(ts['date'], "%Y%m%d") + datetime.timedelta(hours=ts['hour_num'])
            time[ts_idx] = nc.stringtochar(np.array(dt.strftime("%Y-%m-%d %H:%M:%S +000").ljust(64), dtype='S64'))
            
            # Process all species for this time step in parallel
            with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_count()) as executor:
                args = [(ts, spec, upper_spec, all_time_info, static_params, source_locations)
                       for spec, upper_spec in zip(spec_name_str, uppercase_spec_names)]
                for upper_spec, emission_array in executor.map(process_species_time_step, args):
                    emission_vars[upper_spec][ts_idx, :] = emission_array

    print("\nSuccessfully created PALM LOD2 traffic emissions driver with:")
    print(f"- {len(spec_name_str)} species")
    print(f"- {len(time_steps)} time steps")
    print(f"- {len(source_locations)} source locations")
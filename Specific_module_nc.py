import os
import numpy as np
import netCDF4 as nc
import rasterio
from datetime import datetime, timedelta
from Specific_module_config import emis_geotiff_pth, static_pth, static, spec_name_str, active_categories, molar_mass
from Specific_module_main import resample_geotiff_to_palm_grid

def get_time_bands_info(geotiff_path):
    """Get band indices for selected categories."""
    try:
        with rasterio.open(geotiff_path) as src:
            band_descriptions = src.descriptions or [f"Band_{i+1}" for i in range(src.count)]
            selected_indices = [
                i for i, desc in enumerate(band_descriptions)
                if desc and any(desc.startswith(cat) for cat in active_categories)
            ]
            if not selected_indices:
                print(f"Warning: No bands found for {geotiff_path} with categories {active_categories}")
                return []
            if len(selected_indices) != 24:
                print(f"Warning: Found {len(selected_indices)} bands for {geotiff_path}, expected 24")
            for idx in selected_indices[:5]:  # Log sample bands for debugging
                band = src.read(idx + 1)
                valid_vals = band[(~np.isnan(band)) & (band != src.nodatavals[idx])]
                if len(valid_vals) > 0:
                    print(f"  Band {idx + 1} raw values: Min {np.min(valid_vals):.2e}, Max {np.max(valid_vals):.2e}, Unique {len(np.unique(valid_vals))}")
                else:
                    print(f"  Band {idx + 1} raw values: No valid data")
            return selected_indices
    except Exception as e:
        print(f"Error reading GeoTIFF {geotiff_path}: {str(e)}")
        return []

def get_source_locations(static_params, all_time_info):
    """Identify non-zero emission cells across all species and time steps."""
    ny, nx = static_params['ny'], static_params['nx']
    non_zero_threshold = 1e-18
    emission_mask = np.zeros((ny, nx), dtype=bool)

    for spec in spec_name_str:
        gt_path = f"{emis_geotiff_pth}emission_{spec}_temporal.tif"
        if spec not in all_time_info or not all_time_info[spec]:
            continue
        for band_idx in all_time_info[spec]:
            try:
                arr = resample_geotiff_to_palm_grid(gt_path, static_params, band_num=band_idx + 1)
                valid_mask = (~np.isnan(arr)) & (arr != -9999.9) & (np.abs(arr) > non_zero_threshold)
                emission_mask |= valid_mask
                if np.any(valid_mask):
                    valid_vals = arr[valid_mask]
                    print(f"  {spec} band {band_idx + 1} resampled: Min {np.min(valid_vals):.2e}, Max {np.max(valid_vals):.2e}, Unique {len(np.unique(valid_vals))}")
            except Exception as e:
                print(f"Error processing band {band_idx + 1} for {spec}: {str(e)}")
                continue

    i, j = np.nonzero(emission_mask)
    k = np.ones_like(i, dtype=np.int32)  # k=1 as required
    source_locations = [(i[n], j[n], k[n]) for n in range(len(i))]
    
    print(f"Source location statistics:")
    print(f"- Total emission sources: {len(source_locations)}")
    print(f"- Sample source locations: {source_locations[:5] if source_locations else []}")
    return source_locations

def convert_emission_units(arr, spec, dz):
    """Convert emission units from kg/m²/hr to mol/(m³s) for gases or kg/(m³s) for PM."""
    if np.all(arr == 0):
        return arr
    if spec in molar_mass and molar_mass[spec] is not None:
        # Gas phase: kg/m²/hr → mol/m³/s
        converted = arr / (molar_mass[spec] / 1000) / (dz * 3600)
    else:
        # Particulate: kg/m²/hr → kg/m³/s
        converted = arr / (dz * 3600)
    valid_vals = converted[(~np.isnan(converted)) & (converted != -9999.9) & (converted > 0)]
    if len(valid_vals) > 0:
        print(f"Converted {spec}: Min {np.min(valid_vals):.2e}, Max {np.max(valid_vals):.2e}, Unique {len(np.unique(valid_vals))}")
    return converted

def create_chemistry_driver(static_params):
    """Create PALM LOD2 Chemistry Driver for Traffic Emissions."""
    print("\nCreating PALM LOD2 Chemistry Driver for Traffic Emissions")

    # Validate static parameters
    required_keys = ['nx', 'ny', 'dx', 'dy', 'dz', 'origin_x', 'origin_y', 'origin_time']
    if not all(key in static_params for key in required_keys):
        missing = [key for key in required_keys if key not in static_params]
        raise ValueError(f"Missing required static parameters: {missing}")

    # Collect temporal emission data
    all_time_info = {}
    for spec in spec_name_str:
        gt_path = f"{emis_geotiff_pth}emission_{spec}_temporal.tif"
        if os.path.exists(gt_path):
            print(f"Processing {spec}...")
            band_indices = get_time_bands_info(gt_path)
            if band_indices:
                all_time_info[spec] = band_indices
                print(f"  Found {len(band_indices)} bands for {spec}")

    if not all_time_info:
        raise ValueError("No valid emission files found for specified species")

    # Get source locations
    source_locations = get_source_locations(static_params, all_time_info)
    if not source_locations:
        raise ValueError("No valid emission source locations found")

    # Generate timestamps (24 hourly steps from origin_time)
    origin_time = static_params['origin_time'].replace("+00", "+0000")
    try:
        start_dt = datetime.strptime(origin_time, "%Y-%m-%d %H:%M:%S %z")
    except ValueError:
        start_dt = datetime.strptime(origin_time, "%Y-%m-%d %H:%M:%S")
        start_dt = start_dt.replace(tzinfo=pytz.UTC)
    timestamps = [
        (start_dt + timedelta(hours=i)).strftime("%Y-%m-%d %H:%M:%S +000")
        for i in range(24)
    ]

    # Species configuration
    species_mapping = {spec: spec.upper().replace('_', '') for spec in spec_name_str}
    uppercase_spec_names = [species_mapping[spec] for spec in spec_name_str]
    gas_phase_species = {'N2O', 'NOX', 'NMVOC', 'SO2', 'CO', 'NH3', 'CO2', 'CH4', 'NO', 'NO2', 'O3'}

    # Create NetCDF file
    output_path = f"{static_pth}{static}_emis_traffic"
    print(f"Creating output file: {output_path}")
    
    with nc.Dataset(output_path, 'w', format='NETCDF3_CLASSIC') as ds:
        # Global attributes
        ds.setncatts({
            'description': 'LOD2 traffic emissions for PALM model simulation from Downscaled GRETA Emissions',
            'author': 'Sathish Kumar Vaithiyanadhan (sathish.vaithiyanadhan@med.uni-augsburg.de)',
            'Institution': 'Chair of Model-based Environmental Exposure Science (MBEES), Faculty of Medicine, University of Augsburg, Germany.',
            'lod': 2,
            'origin_time': static_params['origin_time'],
            'origin_x': static_params['origin_x'],
            'origin_y': static_params['origin_y'],
            'resolution': static_params['dx'],
            'history': f'Created on {datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S +000")}',
            'source': 'multiple sources; netCDF3',
            'Conventions': 'CF-1.7'
        })

        # Dimensions
        ds.createDimension('ntime', len(timestamps))
        ds.createDimension('nspecies', len(spec_name_str))
        ds.createDimension('nvsrc', len(source_locations))
        ds.createDimension('field_length', 64)

        # Coordinate variables
        timestamp = ds.createVariable('timestamp', 'S1', ('ntime', 'field_length'))
        timestamp.description = 'Time stamps in UTC'
        timestamp[:] = nc.stringtochar(np.array([t.ljust(64) for t in timestamps], dtype='S64'))

        species = ds.createVariable('species', 'S1', ('nspecies', 'field_length'))
        species.description = 'Emission species names'
        species[:] = nc.stringtochar(np.array(uppercase_spec_names, dtype='S64'))

        vsrc_i = ds.createVariable('vsrc_i', 'i4', ('nvsrc',))
        vsrc_i.description = 'i grid indices for volume source location'
        vsrc_i[:] = [loc[0] for loc in source_locations]

        vsrc_j = ds.createVariable('vsrc_j', 'i4', ('nvsrc',))
        vsrc_j.description = 'j grid indices for volume source location'
        vsrc_j[:] = [loc[1] for loc in source_locations]

        vsrc_k = ds.createVariable('vsrc_k', 'i4', ('nvsrc',))
        vsrc_k.description = 'k grid indices for volume source location'
        vsrc_k[:] = [loc[2] for loc in source_locations]

        # Emission variables
        emission_vars = {}
        for spec in uppercase_spec_names:
            var_name = f'vsrc_{spec}'
            units = 'mol/(m^3 s)' if spec in gas_phase_species else 'kg/(m^3 s)'
            emission_vars[spec] = ds.createVariable(var_name, 'f4', ('ntime', 'nvsrc'), fill_value=0.0)
            emission_vars[spec].description = f'Volume source values for {spec}'
            emission_vars[spec].units = units

        # Process emissions
        for ts_idx, ts in enumerate(timestamps):
            print(f"Processing timestep {ts}...")
            for spec, upper_spec in zip(spec_name_str, uppercase_spec_names):
                emission_array = np.zeros(len(source_locations), dtype=np.float32)
                if spec in all_time_info and len(all_time_info[spec]) > ts_idx:
                    band_idx = all_time_info[spec][ts_idx]
                    try:
                        # Read raw GeoTIFF data for debugging
                        with rasterio.open(f"{emis_geotiff_pth}emission_{spec}_temporal.tif") as src:
                            raw_data = src.read(band_idx + 1)
                            valid_raw = raw_data[(~np.isnan(raw_data)) & (raw_data != src.nodatavals[band_idx])]
                            if len(valid_raw) > 0:
                                print(f"  {spec} band {band_idx + 1} raw: Min {np.min(valid_raw):.2e}, Max {np.max(valid_raw):.2e}, Unique {len(np.unique(valid_raw))}")
                            else:
                                print(f"  {spec} band {band_idx + 1} raw: No valid data")
                        
                        arr = resample_geotiff_to_palm_grid(
                            f"{emis_geotiff_pth}emission_{spec}_temporal.tif",
                            static_params,
                            band_num=band_idx + 1
                        )
                        # Log resampled data
                        valid_resampled = arr[(~np.isnan(arr)) & (arr != -9999.9)]
                        if len(valid_resampled) > 0:
                            print(f"  {spec} band {band_idx + 1} resampled: Min {np.min(valid_resampled):.2e}, Max {np.max(valid_resampled):.2e}, Unique {len(np.unique(valid_resampled))}")
                        
                        arr = convert_emission_units(arr, spec, static_params['dz'])
                        valid_vals = arr[(~np.isnan(arr)) & (arr != -9999.9) & (arr > 0)]
                        if len(valid_vals) > 0:
                            print(f"  {spec} timestep {ts_idx + 1}: Min {np.min(valid_vals):.2e}, Max {np.max(valid_vals):.2e}, Unique {len(np.unique(valid_vals))}")
                        else:
                            print(f"  {spec} timestep {ts_idx + 1}: No valid values")
                        
                        for src_idx, (i, j, k) in enumerate(source_locations):
                            if 0 <= j < static_params['ny'] and 0 <= i < static_params['nx']:
                                emission_array[src_idx] = arr[j, i] if not np.isnan(arr[j, i]) else 0.0
                    except Exception as e:
                        print(f"Error processing {spec} band {band_idx + 1}: {str(e)}")
                        emission_array.fill(0.0)
                emission_vars[upper_spec][ts_idx, :] = emission_array

    print("\nSuccessfully created emissions driver with:")
    print(f"- {len(spec_name_str)} species")
    print(f"- {len(timestamps)} time steps")
    print(f"- {len(source_locations)} source locations")
    print(f"Output file: {output_path}")

    # Final verification
    print("\nFinal verification:")
    with nc.Dataset(output_path) as ds:
        for var in ds.variables:
            if var.startswith('vsrc_') and var not in ['vsrc_i', 'vsrc_j', 'vsrc_k']:
                data = ds[var][:]
                valid_data = data[data > 0]
                print(f"\nVariable {var}:")
                print(f"  Shape: {data.shape}")
                print(f"  Non-zero values: {len(valid_data)}")
                if len(valid_data) > 0:
                    print(f"  Min: {np.min(valid_data):.2e}")
                    print(f"  Max: {np.max(valid_data):.2e}")
                    print(f"  Mean: {np.mean(valid_data):.2e}")
                    print(f"  Unique values: {len(np.unique(valid_data))}")
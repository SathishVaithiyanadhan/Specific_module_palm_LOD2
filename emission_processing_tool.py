

## working
'''import rasterio
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import xarray as xr
from datetime import datetime, timedelta
import pytz
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds
import os
import netCDF4 as nc
from emission_config import geotiff_dir, static_path, output_dir, selected_band_prefix, active_categories, spec_name_str, molar_mass, layer_height, non_zero_threshold, model_name, emission_mode, field_length

# Define species name mapping
species_mapping = {
    'n2o': 'N2O', 'nox': 'NOX', 'nmvoc': 'NMVOC', 'so2': 'SO2', 'co': 'CO',
    'pm10': 'PM10', 'pm2_5': 'PM25', 'nh3': 'NH3', 'pb': 'PB', 'cd': 'CD',
    'hg': 'HG', 'as': 'AS', 'ni': 'NI', 'bc': 'BC', 'co2': 'CO2', 'ch4': 'CH4',
    'no': 'NO', 'no2': 'NO2', 'ec': 'EC', 'na': 'NA', 'so4': 'SO4', 'oc': 'OC',
    'othmin': 'OTHMIN', 'o3': 'O3'
}

# Step 1: Read static NetCDF data
try:
    with xr.open_dataset(static_path) as ds_static:
        origin_x = ds_static.attrs["origin_x"]
        origin_y = ds_static.attrs["origin_y"]
        origin_time = ds_static.attrs["origin_time"]
        x = ds_static["x"].values
        y = ds_static["y"].values
        z = ds_static["z"].values
        nx, ny, nz = len(x), len(y), len(z)
        dx = ds_static.attrs.get("dx", 3.0)  # Get resolution from attributes or default to 3.0
except Exception as e:
    raise Exception(f"Failed to read static NetCDF: {e}")

# Compute geographic extent
extent = [origin_x, origin_x + nx * dx, origin_y, origin_y + ny * dx]
print(f"Static data extent: left={extent[0]}, right={extent[1]}, bottom={extent[2]}, top={extent[3]}")
print(f"Grid size: nx={nx}, ny={ny}, nz={nz}")
print(f"Origin time: {origin_time}")

# Step 2: Read GeoTIFF and process all categories and species
vsrc_data = {}
n_bands = 0
selected_categories = [c for c in active_categories if c.startswith(selected_band_prefix) or selected_band_prefix == "SumAllSectors"]
selected_species = [species_mapping[spec] for spec in spec_name_str if spec in species_mapping]
unmapped_species = [spec for spec in spec_name_str if spec not in species_mapping]
if unmapped_species:
    print(f"Warning: The following species are not in species_mapping and will be skipped: {unmapped_species}")
molar_mass_upper = {species_mapping.get(k, k.upper()): v for k, v in molar_mass.items()}

for species in selected_species:
    input_species = next((k for k, v in species_mapping.items() if v == species), species.lower())
    emissions_sum = None
    unit = "mol/(m^3 s)" if molar_mass_upper.get(species) is not None else "kg/(m^3 s)"
    
    for category in selected_categories:
        geotiff_path = f"{geotiff_dir}emission_{input_species}_temporal.tif"
        if not os.path.exists(geotiff_path):
            print(f"GeoTIFF file not found: {geotiff_path}")
            continue
            
        try:
            with rasterio.open(geotiff_path) as src:
                print(f"GeoTIFF {geotiff_path}: Bounds {src.bounds}, Resolution {src.res}")
                band_descriptions = src.descriptions or [f"Band_{i+1}" for i in range(src.count)]
                selected_band_indices = [
                    i for i, desc in enumerate(band_descriptions)
                    if desc and desc.startswith(category)
                ]
                if not selected_band_indices:
                    continue
                    
                emissions_kg_m2_hr = src.read([i + 1 for i in selected_band_indices])
                n_bands, rows, cols = emissions_kg_m2_hr.shape
                src_transform = src.transform
                src_crs = src.crs
                # Log raw data
                valid_raw = emissions_kg_m2_hr[(~np.isnan(emissions_kg_m2_hr)) & (emissions_kg_m2_hr != src.nodatavals[0])]
                if len(valid_raw) > 0:
                    print(f"{species} ({category}) raw: Min {valid_raw.min():.2e}, Max {valid_raw.max():.2e}")

                # Resample to static grid
                target_transform = from_bounds(extent[0], extent[2], extent[1], extent[3], nx, ny)
                emissions_resampled = np.zeros((n_bands, ny, nx), dtype=np.float32)
                for band in range(n_bands):
                    reproject(
                        source=emissions_kg_m2_hr[band],
                        destination=emissions_resampled[band],
                        src_transform=src_transform,
                        src_crs=src_crs,
                        dst_transform=target_transform,
                        dst_crs=src_crs,
                        resampling=Resampling.bilinear
                    )
                    valid_resampled = emissions_resampled[band][~np.isnan(emissions_resampled[band])]
                    if len(valid_resampled) > 0:
                        print(f"{species} ({category}) Band {band + 1} resampled: Min {valid_resampled.min():.2e}, Max {valid_resampled.max():.2e}")

                # Unit conversion
                emissions_kg_m3_s = emissions_resampled / (layer_height * 3600)
                if molar_mass_upper.get(species) is not None:
                    emissions_converted = emissions_kg_m3_s * 1000 / molar_mass_upper[species]
                else:
                    emissions_converted = emissions_kg_m3_s
                valid_converted = emissions_converted[~np.isnan(emissions_converted)]
                if len(valid_converted) > 0:
                    print(f"{species} ({category}) converted: Min {valid_converted.min():.2e}, Max {valid_converted.max():.2e}")

                if emissions_sum is None:
                    emissions_sum = emissions_converted
                else:
                    emissions_sum += emissions_converted

        except Exception as e:
            print(f"Error processing {species} ({category}): {e}")
            continue
    
    if emissions_sum is not None:
        vsrc_data[species] = (emissions_sum, unit)

# Check for valid data
if not vsrc_data:
    raise ValueError("No valid emission data found for any species or category.")

# Step 3: Compute volumetric emission sources and cell coordinates
has_emission = None
for species, (data, _) in vsrc_data.items():
    has_emission_temp = np.any(np.abs(data) > non_zero_threshold, axis=0)
    has_emission = has_emission_temp if has_emission is None else np.logical_or(has_emission, has_emission_temp)

i, j = np.where(has_emission)
k = np.ones_like(i, dtype=np.int32)  # k=1 (surface layer)
vsrc_i, vsrc_j, vsrc_k = i, j, k
print(f"Number of source locations: {len(vsrc_i)}")
print(f"Sample source locations: {list(zip(vsrc_i[:5], vsrc_j[:5], vsrc_k[:5]))}")
vsrc_species = {f"vsrc_{species}": data[:, vsrc_i, vsrc_j] for species, (data, _) in vsrc_data.items()}

# Step 4: Generate timestamps
origin_time_fixed = origin_time.replace("+00", "+0000")
start_datetime = datetime.strptime(origin_time_fixed, "%Y-%m-%d %H:%M:%S %z")
timestamps = [
    (start_datetime + timedelta(hours=i)).strftime("%Y-%m-%d %H:%M:%S +000")
    for i in range(n_bands)
]

# Step 5: Create NetCDF file following PALM specifications
output_path = f"{output_dir}{model_name}_emis_traffic"  # No .nc extension
print(f"Creating output file: {output_path}")
with nc.Dataset(output_path, 'w', format='NETCDF3_CLASSIC') as ds:
    # Create required dimensions
    ds.createDimension('ntime', n_bands)
    ds.createDimension('nspecies', len(vsrc_data))
    ds.createDimension('nvsrc', len(vsrc_i))
    ds.createDimension('field_length', field_length)

    # Create variables exactly as specified in PALM documentation
    timestamp_var = ds.createVariable('timestamp', 'S1', ('ntime', 'field_length'))
    timestamp_var.description = "Time stamps in UTC"
    timestamp_var[:] = nc.stringtochar(np.array([t.ljust(field_length) for t in timestamps], dtype='S64'))

    species_var = ds.createVariable('species', 'S1', ('nspecies', 'field_length'))
    species_var.description = "Emission species names"
    species_list = list(vsrc_data.keys())
    species_var[:] = nc.stringtochar(np.array([s.ljust(field_length) for s in species_list], dtype='S64'))

    vsrc_i_var = ds.createVariable('vsrc_i', 'i4', ('nvsrc',))
    vsrc_i_var.description = "i grid indices for volume source location"
    vsrc_i_var[:] = vsrc_i

    vsrc_j_var = ds.createVariable('vsrc_j', 'i4', ('nvsrc',))
    vsrc_j_var.description = "j grid indices for volume source location"
    vsrc_j_var[:] = vsrc_j

    vsrc_k_var = ds.createVariable('vsrc_k', 'i4', ('nvsrc',))
    vsrc_k_var.description = "k grid indices for volume source location"
    vsrc_k_var[:] = vsrc_k

    for species in species_list:
        var_name = f"vsrc_{species}"
        unit = vsrc_data[species][1]
        var = ds.createVariable(var_name, 'f4', ('ntime', 'nvsrc'), fill_value=0.0)
        var.description = f"Volume source values for {species}"
        var.units = unit
        var[:] = vsrc_species[var_name]

    ds.setncatts({
        'description': 'LOD2 traffic emissions for PALM model simulation from Downscaled GRETA Emissions',
        'author': 'Sathish Kumar Vaithiyanadhan (sathish.vaithiyanadhan@med.uni-augsburg.de)',
        'Institution': 'Chair of Model-based Environmental Exposure Science (MBEES), Faculty of Medicine, University of Augsburg, Germany.',
        'lod': 2,
        'origin_time': origin_time,
        'origin_x': origin_x,
        'origin_y': origin_y,
        'resolution': dx,
        'history': f'Created on {datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S +000")}',
        'source': 'multiple sources; netCDF3',
        'Conventions': 'CF-1.7'
    })

print(f"Successfully created PALM-compatible NetCDF file: {output_path}")

# Step 6: Final verification
with nc.Dataset(output_path) as ds:
    for var in ds.variables:
        if var.startswith('vsrc_') and var not in ['vsrc_i', 'vsrc_j', 'vsrc_k']:
            data = ds[var][:]
            valid_data = data[data > 0]
            print(f"\nVariable {var}:")
            print(f"  Shape: {data.shape}")
            print(f"  Non-zero values: {len(valid_data)}/{data.size} ({len(valid_data)/data.size*100:.2f}%)")
            if len(valid_data) > 0:
                print(f"  Min: {np.min(valid_data):.2e}")
                print(f"  Max: {np.max(valid_data):.2e}")
                print(f"  Mean: {np.mean(valid_data):.2e}")
                max_idx = np.unravel_index(np.argmax(data, axis=None), data.shape)
                print(f"  Peak emission at: Timestep {max_idx[0] + 1} - {timestamps[max_idx[0]]}")

# Step 7: Visualization (optional)
if src_crs and src_crs.to_epsg() == 25832:
    projection = ccrs.epsg(25832)
else:
    projection = ccrs.PlateCarree()

for species, (data, unit) in vsrc_data.items():
    for t in range(min(3, n_bands)):
        fig = plt.figure(figsize=(10, 8))
        ax = plt.axes(projection=projection)
        emissions = np.ma.masked_invalid(data[t, :, :])
        im = ax.imshow(emissions, origin="upper", extent=extent, cmap="viridis", transform=projection)
        plt.colorbar(im, label=f"{species} Emissions ({unit})")
        plt.title(f"{species} Emissions at {timestamps[t]}")
        #plt.show()
        plt.close(fig)'''

###
## flipped
import rasterio
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import xarray as xr
from datetime import datetime, timedelta
import pytz
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds
import os
import netCDF4 as nc
from emission_config import geotiff_dir, static_path, output_dir, selected_band_prefix, active_categories, spec_name_str, molar_mass, layer_height, non_zero_threshold, model_name, emission_mode, field_length

# Define species name mapping
species_mapping = {
    'n2o': 'N2O', 'nox': 'NOX', 'nmvoc': 'NMVOC', 'so2': 'SO2', 'co': 'CO',
    'pm10': 'PM10', 'pm2_5': 'PM25', 'nh3': 'NH3', 'pb': 'PB', 'cd': 'CD',
    'hg': 'HG', 'as': 'AS', 'ni': 'NI', 'bc': 'BC', 'co2': 'CO2', 'ch4': 'CH4',
    'no': 'NO', 'no2': 'NO2', 'ec': 'EC', 'na': 'NA', 'so4': 'SO4', 'oc': 'OC',
    'othmin': 'OTHMIN', 'o3': 'O3'
}

# Step 1: Read static NetCDF data
try:
    with xr.open_dataset(static_path) as ds_static:
        origin_x = ds_static.attrs["origin_x"]
        origin_y = ds_static.attrs["origin_y"]
        origin_time = ds_static.attrs["origin_time"]
        x = ds_static["x"].values
        y = ds_static["y"].values
        z = ds_static["z"].values
        nx, ny, nz = len(x), len(y), len(z)
        dx = ds_static.attrs.get("dx", 3.0)
except Exception as e:
    raise Exception(f"Failed to read static NetCDF: {e}")

# Verify and adjust grid orientation
if not np.all(np.diff(x) > 0):
    raise ValueError("Static NetCDF x-coordinates are not in ascending order (required for east-right orientation).")
if not np.all(np.diff(y) > 0):
    print("Warning: Static NetCDF y-coordinates are not in ascending order. Flipping y-axis to ensure north-up orientation.")
    y = y[::-1]
    origin_y = y[0]
    flip_y = True
else:
    flip_y = False

# Compute geographic extent (adjusted for transposition)
extent = [origin_x, origin_x + nx * dx, origin_y, origin_y + ny * dx]
print(f"Static data extent: left={extent[0]}, right={extent[1]}, bottom={extent[2]}, top={extent[3]}")
print(f"Grid size: nx={nx}, ny={ny}, nz={nz}")
print(f"Origin time: {origin_time}")

# Step 2: Read GeoTIFF and process all categories and species
vsrc_data = {}
n_bands = 0
selected_categories = [c for c in active_categories if c.startswith(selected_band_prefix) or selected_band_prefix == "SumAllSectors"]
selected_species = [species_mapping[spec] for spec in spec_name_str if spec in species_mapping]
unmapped_species = [spec for spec in spec_name_str if spec not in species_mapping]
if unmapped_species:
    print(f"Warning: The following species are not in species_mapping and will be skipped: {unmapped_species}")
molar_mass_upper = {species_mapping.get(k, k.upper()): v for k, v in molar_mass.items()}

# Initialize src_crs to None; will be set from the first valid GeoTIFF
src_crs = None

for species in selected_species:
    input_species = next((k for k, v in species_mapping.items() if v == species), species.lower())
    emissions_sum = None
    unit = "mol/(m^3 s)" if molar_mass_upper.get(species) is not None else "kg/(m^3 s)"
    
    for category in selected_categories:
        geotiff_path = f"{geotiff_dir}emission_{input_species}_temporal.tif"
        if not os.path.exists(geotiff_path):
            print(f"GeoTIFF file not found: {geotiff_path}")
            continue
            
        try:
            with rasterio.open(geotiff_path) as src:
                print(f"GeoTIFF {geotiff_path}: Bounds {src.bounds}, Resolution {src.res}")
                if src_crs is None:
                    src_crs = src.crs
                elif src.crs != src_crs:
                    print(f"Warning: CRS mismatch in {geotiff_path}. Expected {src_crs}, got {src.crs}")
                
                band_descriptions = src.descriptions or [f"Band_{i+1}" for i in range(src.count)]
                selected_band_indices = [
                    i for i, desc in enumerate(band_descriptions)
                    if desc and desc.startswith(category)
                ]
                if not selected_band_indices:
                    continue
                emissions_kg_m2_hr = src.read([i + 1 for i in selected_band_indices])
                n_bands, rows, cols = emissions_kg_m2_hr.shape
                src_transform = src.transform
                valid_raw = emissions_kg_m2_hr[(~np.isnan(emissions_kg_m2_hr)) & (emissions_kg_m2_hr != src.nodatavals[0])]
                if len(valid_raw) > 0:
                    print(f"{species} ({category}) raw: Min {valid_raw.min():.2e}, Max {valid_raw.max():.2e}")

                # Resample to static grid
                target_transform = from_bounds(extent[0], extent[2], extent[1], extent[3], nx, ny)
                emissions_resampled = np.zeros((n_bands, ny, nx), dtype=np.float32)
                for band in range(n_bands):
                    reproject(
                        source=emissions_kg_m2_hr[band],
                        destination=emissions_resampled[band],
                        src_transform=src_transform,
                        src_crs=src_crs,
                        dst_transform=target_transform,
                        dst_crs=src_crs,
                        resampling=Resampling.bilinear
                    )
                    valid_resampled = emissions_resampled[band][~np.isnan(emissions_resampled[band])]
                    if len(valid_resampled) > 0:
                        print(f"{species} ({category}) Band {band + 1} resampled: Min {valid_resampled.min():.2e}, Max {valid_resampled.max():.2e}")

                # Flip y-axis if static grid was flipped
                if flip_y:
                    emissions_resampled = emissions_resampled[:, ::-1, :]

                # Transpose to rotate 90 degrees left (swap y and x axes)
                emissions_resampled = emissions_resampled.transpose(0, 2, 1)  # Shape becomes (n_bands, nx, ny)

                # Unit conversion
                emissions_kg_m3_s = emissions_resampled / (layer_height * 3600)
                if molar_mass_upper.get(species) is not None:
                    emissions_converted = emissions_kg_m3_s * 1000 / molar_mass_upper[species]
                else:
                    emissions_converted = emissions_kg_m3_s
                valid_converted = emissions_converted[~np.isnan(emissions_converted)]
                if len(valid_converted) > 0:
                    print(f"{species} ({category}) converted: Min {valid_converted.min():.2e}, Max {valid_converted.max():.2e}")

                if emissions_sum is None:
                    emissions_sum = emissions_converted
                else:
                    emissions_sum += emissions_converted

        except Exception as e:
            print(f"Error processing {species} ({category}): {e}")
            continue
    
    if emissions_sum is not None:
        vsrc_data[species] = (emissions_sum, unit)

# Check for valid data
if not vsrc_data:
    raise ValueError("No valid emission data found for any species or category.")

# Step 3: Compute volumetric emission sources and cell coordinates
has_emission = None
for species, (data, _) in vsrc_data.items():
    has_emission_temp = np.any(np.abs(data) > non_zero_threshold, axis=0)
    has_emission = has_emission_temp if has_emission is None else np.logical_or(has_emission, has_emission_temp)

# Adjust indices for transposed orientation
j, i = np.where(has_emission)  # Swap i and j to match transposed array
if flip_y:
    i = nx - 1 - i  # Adjust i for north-up in transposed grid
# Swap i and j to correct orientation (vsrc_i as x, vsrc_j as y)
vsrc_i, vsrc_j, vsrc_k = j, i, np.ones_like(i, dtype=np.int32)  # Reversed to match PALM's expected layout
print(f"Number of source locations: {len(vsrc_i)}")
print(f"Sample source locations: {list(zip(vsrc_i[:5], vsrc_j[:5], vsrc_k[:5]))}")
vsrc_species = {f"vsrc_{species}": data[:, vsrc_i, vsrc_j] for species, (data, _) in vsrc_data.items()}

# Step 4: Generate timestamps
origin_time_fixed = origin_time.replace("+00", "+0000")
start_datetime = datetime.strptime(origin_time_fixed, "%Y-%m-%d %H:%M:%S %z")
timestamps = [
    (start_datetime + timedelta(hours=i)).strftime("%Y-%m-%d %H:%M:%S +000")
    for i in range(n_bands)
]

# Step 5: Create NetCDF file following PALM specifications
output_path = f"{output_dir}{model_name}_emis_traffic"
print(f"Creating output file: {output_path}")
with nc.Dataset(output_path, 'w', format='NETCDF3_CLASSIC') as ds:
    # Create required dimensions (swap nx and ny)
    ds.createDimension('ntime', n_bands)
    ds.createDimension('nspecies', len(vsrc_data))
    ds.createDimension('nvsrc', len(vsrc_i))
    ds.createDimension('field_length', field_length)
    ds.createDimension('x', nx)
    ds.createDimension('y', ny)

    # Create grid mapping variable for CRS
    if src_crs is not None:
        grid_mapping = ds.createVariable('crs', 'i4')
        grid_mapping.grid_mapping_name = "latitude_longitude" if src_crs.is_geographic else "transverse_mercator"
        if src_crs.to_epsg():
            grid_mapping.epsg_code = f"EPSG:{src_crs.to_epsg()}"
        grid_mapping.crs_wkt = src_crs.to_wkt()
        grid_mapping.description = "Coordinate Reference System from input GeoTIFF"

    # Create coordinate variables (swap x and y roles due to transposition)
    x_var = ds.createVariable('x', 'f4', ('x',))  # Vertical axis (original y, northing)
    x_var.units = "meter"
    x_var.description = "Northing coordinates in EPSG:25832"
    x_var[:] = y  # Use y for northing

    y_var = ds.createVariable('y', 'f4', ('y',))  # Horizontal axis (original x, easting)
    y_var.units = "meter"
    y_var.description = "Easting coordinates in EPSG:25832"
    y_var[:] = x  # Use x for easting

    # Create variables
    timestamp_var = ds.createVariable('timestamp', 'S1', ('ntime', 'field_length'))
    timestamp_var.description = "Time stamps in UTC"
    timestamp_var[:] = nc.stringtochar(np.array([t.ljust(field_length) for t in timestamps], dtype='S64'))

    species_var = ds.createVariable('species', 'S1', ('nspecies', 'field_length'))
    species_var.description = "Emission species names"
    species_list = list(vsrc_data.keys())
    species_var[:] = nc.stringtochar(np.array([s.ljust(field_length) for s in species_list], dtype='S64'))

    vsrc_i_var = ds.createVariable('vsrc_i', 'i4', ('nvsrc',))
    vsrc_i_var.description = "i grid indices for volume source location (x-axis, north-up)"
    vsrc_i_var.grid_mapping = "crs" if src_crs is not None else ""
    vsrc_i_var.coordinates = "x y"
    vsrc_i_var[:] = vsrc_i

    vsrc_j_var = ds.createVariable('vsrc_j', 'i4', ('nvsrc',))
    vsrc_j_var.description = "j grid indices for volume source location (y-axis, east-right)"
    vsrc_j_var.grid_mapping = "crs" if src_crs is not None else ""
    vsrc_j_var.coordinates = "x y"
    vsrc_j_var[:] = vsrc_j

    vsrc_k_var = ds.createVariable('vsrc_k', 'i4', ('nvsrc',))
    vsrc_k_var.description = "k grid indices for volume source location"
    vsrc_k_var.grid_mapping = "crs" if src_crs is not None else ""
    vsrc_k_var.coordinates = "x y"
    vsrc_k_var[:] = vsrc_k

    for species in species_list:
        var_name = f"vsrc_{species}"
        unit = vsrc_data[species][1]
        var = ds.createVariable(var_name, 'f4', ('ntime', 'nvsrc'), fill_value=0.0)
        var.description = f"Volume source values for {species}"
        var.units = unit
        var.grid_mapping = "crs" if src_crs is not None else ""
        var.coordinates = "x y"
        var[:] = vsrc_species[var_name]

    ds.setncatts({
        'description': 'LOD2 traffic emissions for PALM model simulation from Downscaled GRETA Emissions',
        'author': 'Sathish Kumar Vaithiyanadhan (sathish.vaithiyanadhan@med.uni-augsburg.de)',
        'Institution': 'Chair of Model-based Environmental Exposure Science (MBEES), Faculty of Medicine, University of Augsburg, Germany.',
        'lod': 2,
        'origin_time': origin_time,
        'origin_x': origin_x,
        'origin_y': origin_y,
        'resolution': dx,
        'history': f'Created on {datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S +000")}',
        'source': 'multiple sources; netCDF3',
        'Conventions': 'CF-1.7'
    })

print(f"Successfully created PALM-compatible NetCDF file: {output_path}")

# Step 6: Final verification
with nc.Dataset(output_path) as ds:
    for var in ds.variables:
        if var.startswith('vsrc_') and var not in ['vsrc_i', 'vsrc_j', 'vsrc_k']:
            data = ds[var][:]
            valid_data = data[data > 0]
            print(f"\nVariable {var}:")
            print(f"  Shape: {data.shape}")
            print(f"  Non-zero values: {len(valid_data)}/{data.size} ({len(valid_data)/data.size*100:.2f}%)")
            if len(valid_data) > 0:
                print(f"  Min: {np.min(valid_data):.2e}")
                print(f"  Max: {np.max(valid_data):.2e}")
                print(f"  Mean: {np.mean(valid_data):.2e}")
                max_idx = np.unravel_index(np.argmax(data, axis=None), data.shape)
                print(f"  Peak emission at: Timestep {max_idx[0] + 1} - {timestamps[max_idx[0]]}")

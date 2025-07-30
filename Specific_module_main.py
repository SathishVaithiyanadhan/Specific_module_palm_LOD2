import netCDF4 as nc
import numpy as np
from osgeo import gdal, osr
from Specific_module_config import transformer_to_utm, transformer_to_wgs, config_proj, static_pth, static

def extract_static_parameters(static_file):
    """Extract domain parameters from static driver file dynamically"""
    with nc.Dataset(static_file, "r") as ncs:
        params = {}
        # Grid dimensions
        params['nx'] = ncs.dimensions['x'].size if 'x' in ncs.dimensions else 1
        params['ny'] = ncs.dimensions['y'].size if 'y' in ncs.dimensions else 1
        params['nz'] = ncs.dimensions['z'].size if 'z' in ncs.dimensions else 1
        # Grid resolution (use res_orig from buildings_2d)
        if 'buildings_2d' not in ncs.variables or 'res_orig' not in ncs.variables['buildings_2d'].ncattrs():
            raise ValueError("Static file missing 'buildings_2d' variable or 'res_orig' attribute")
        params['dx'] = ncs.variables['buildings_2d'].getncattr('res_orig')
        params['dy'] = params['dx']  # Assume uniform horizontal resolution
        # Validate res_orig consistency across variables
        for var in ['pavement_type', 'surface_fraction', 'building_type']:
            if var in ncs.variables and 'res_orig' in ncs.variables[var].ncattrs():
                if not np.isclose(ncs.variables[var].getncattr('res_orig'), params['dx'], rtol=1e-5):
                    raise ValueError(f"Inconsistent res_orig in {var}: expected {params['dx']}, found {ncs.variables[var].getncattr('res_orig')}")
        # Origin coordinates (center of domain)
        params['origin_x'] = ncs.getncattr('origin_x') if 'origin_x' in ncs.ncattrs() else 0.0
        params['origin_y'] = ncs.getncattr('origin_y') if 'origin_y' in ncs.ncattrs() else 0.0
        params['origin_lat'] = ncs.getncattr('origin_lat') if 'origin_lat' in ncs.ncattrs() else 0.0
        params['origin_lon'] = ncs.getncattr('origin_lon') if 'origin_lon' in ncs.ncattrs() else 0.0
        params['origin_time'] = ncs.getncattr('origin_time') if 'origin_time' in ncs.ncattrs() else '1970-01-01 00:00:00'
        # Calculate west and south boundaries (lower-left corner)
        params['west'] = params['origin_x'] - (params['nx'] * params['dx']) / 2
        params['south'] = params['origin_y'] - (params['ny'] * params['dy']) / 2
        # Vertical coordinates
        if 'z' not in ncs.variables:
            raise ValueError("Static file missing 'z' variable; cannot compute dz")
        z_coords = ncs.variables['z'][:]
        if len(z_coords) < 2:
            raise ValueError("Static file 'z' variable has fewer than 2 values; cannot compute dz")
        params['z_coords'] = z_coords
        params['dz'] = z_coords[1] - z_coords[0]  # dz for the lowest layer (k=0)
        params['z_origin'] = z_coords[0]
        # Check for non-uniform grid
        dz_array = np.diff(z_coords)
        print(f"Vertical coordinates (z): {z_coords}")
        if not np.allclose(dz_array, dz_array[0], rtol=1e-5):
            print(f"Warning: Non-uniform vertical grid detected (dz varies: {dz_array}). Using dz={params['dz']}m for surface emissions.")
        else:
            print(f"Uniform vertical grid confirmed: dz={params['dz']}m")
        # Building heights
        params['building_height'] = ncs.variables['buildings_2d'][:] if 'buildings_2d' in ncs.variables else np.zeros((params['ny'], params['nx']))
    return params

def read_geotiff_band(geotiff_path, band_num, static_params):
    """Read and resample a GeoTIFF band to match PALM grid"""
    ds = gdal.Open(geotiff_path)
    if not ds:
        raise ValueError(f"Cannot open GeoTIFF: {geotiff_path}")
    band = ds.GetRasterBand(band_num)
    arr = band.ReadAsArray()
    # Replace no-data values
    no_data_value = band.GetNoDataValue() if band.GetNoDataValue() is not None else -9999.9
    arr = np.where(arr == no_data_value, np.nan, arr)
    # Resample to PALM grid
    resampled_arr = resample_geotiff(geotiff_path, static_params)
    ds = None
    return resampled_arr

def resample_geotiff(geotiff_path, static_params):
    """Resample GeoTIFF to match PALM grid resolution and extent"""
    ds = gdal.Open(geotiff_path)
    if not ds:
        raise ValueError(f"Cannot open GeoTIFF: {geotiff_path}")
    
    # Get GeoTIFF metadata
    gt = ds.GetGeoTransform()
    src_cols, src_rows = ds.RasterXSize, ds.RasterYSize
    src_res = gt[1]  # Assuming square pixels (xres = yres)
    
    # Target grid parameters (PALM grid)
    nx, ny = static_params['nx'], static_params['ny']
    dx, dy = static_params['dx'], static_params['dy']
    west, south = static_params['west'], static_params['south']
    east = west + nx * dx
    north = south + ny * dy
    
    # Create in-memory dataset for resampling
    driver = gdal.GetDriverByName('MEM')
    target_ds = driver.Create('', nx, ny, 1, gdal.GDT_Float32)
    
    # Set target geotransform
    target_gt = (west, dx, 0, north, 0, -dy)
    target_ds.SetGeoTransform(target_gt)
    
    # Set projection (assume same as input GeoTIFF, EPSG:25832)
    target_ds.SetProjection(ds.GetProjection())
    
    # Resample
    gdal.ReprojectImage(
        ds,
        target_ds,
        ds.GetProjection(),
        target_ds.GetProjection(),
        gdal.GRA_Bilinear
    )
    
    # Read resampled data
    resampled_arr = target_ds.GetRasterBand(1).ReadAsArray()
    target_ds = None
    ds = None
    return resampled_arr

if __name__ == "__main__":
    from Specific_module_nc import create_chemistry_driver
    print("\nExtracting static driver parameters")
    static_params = extract_static_parameters(f"{static_pth}{static}_static")
    print("\nStatic Driver Configuration:")
    print(f"Center (WGS84): {static_params['origin_lat']}°N, {static_params['origin_lon']}°E")
    print(f"Center (UTM): {static_params['origin_x']}m, {static_params['origin_y']}m")
    print(f"Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")
    print(f"Resolution: {static_params['dx']}m x {static_params['dy']}m x {static_params['dz']}m")
    print(f"Z levels: {static_params['nz']} from {static_params['z_origin']}m")
    #print("\nCreating PALM LOD2 Chemistry Driver for Traffic Emissions")
    create_chemistry_driver(static_params)
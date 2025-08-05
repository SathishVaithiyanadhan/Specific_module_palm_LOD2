import netCDF4 as nc
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_bounds
from Specific_module_config import static_pth, static

def extract_static_parameters(static_file):
    """Extract domain parameters from static driver file."""
    with nc.Dataset(static_file, "r") as ncs:
        params = {}
        # Grid dimensions
        params['nx'] = ncs.dimensions['x'].size if 'x' in ncs.dimensions else 1
        params['ny'] = ncs.dimensions['y'].size if 'y' in ncs.dimensions else 1
        params['nz'] = ncs.dimensions['z'].size if 'z' in ncs.dimensions else 1
        # Grid resolution
        if 'buildings_2d' not in ncs.variables or 'res_orig' not in ncs.variables['buildings_2d'].ncattrs():
            raise ValueError("Static file missing 'buildings_2d' variable or 'res_orig' attribute")
        params['dx'] = ncs.variables['buildings_2d'].getncattr('res_orig')
        params['dy'] = params['dx']  # Assume uniform horizontal resolution
        # Origin coordinates
        params['origin_x'] = ncs.getncattr('origin_x') if 'origin_x' in ncs.ncattrs() else 0.0
        params['origin_y'] = ncs.getncattr('origin_y') if 'origin_y' in ncs.ncattrs() else 0.0
        params['origin_time'] = ncs.getncattr('origin_time') if 'origin_time' in ncs.ncattrs() else '1970-01-01 00:00:00'
        # Vertical coordinates
        if 'z' not in ncs.variables:
            raise ValueError("Static file missing 'z' variable")
        params['z'] = ncs.variables['z'][:]
        if len(params['z']) < 2:
            raise ValueError("Static file 'z' variable has fewer than 2 values")
        params['dz'] = params['z'][1] - params['z'][0]  # dz for the lowest layer
        return params

def resample_geotiff_to_palm_grid(geotiff_path, static_params, band_num=1):
    """Resample GeoTIFF to PALM grid resolution and extent."""
    try:
        with rasterio.open(geotiff_path) as src:
            # Get GeoTIFF metadata
            band = src.read(band_num)
            no_data_value = src.nodatavals[band_num-1] if src.nodatavals[band_num-1] is not None else -9999.9
            band = np.where(band == no_data_value, np.nan, band)
            src_transform = src.transform
            src_crs = src.crs
            src_rows, src_cols = band.shape

            # Target grid parameters (PALM grid)
            nx, ny = static_params['nx'], static_params['ny']
            dx, dy = static_params['dx'], static_params['dy']
            left = static_params['origin_x'] - (nx * dx) / 2
            right = left + nx * dx
            bottom = static_params['origin_y'] - (ny * dy) / 2
            top = bottom + ny * dy
            target_transform = from_bounds(left, bottom, right, top, nx, ny)

            # Resample
            resampled = np.zeros((ny, nx), dtype=np.float32)
            reproject(
                source=band,
                destination=resampled,
                src_transform=src_transform,
                src_crs=src_crs,
                dst_transform=target_transform,
                dst_crs=src_crs,
                resampling=Resampling.bilinear
            )
            return resampled
    except Exception as e:
        print(f"Error resampling GeoTIFF band {band_num}: {str(e)}")
        return np.full((static_params['ny'], static_params['nx']), np.nan)

if __name__ == "__main__":
    from Specific_module_nc import create_chemistry_driver
    print("\nExtracting static driver parameters")
    static_params = extract_static_parameters(f"{static_pth}{static}_static")
    print("\nStatic Driver Configuration:")
    print(f"Center (UTM): {static_params['origin_x']}m, {static_params['origin_y']}m")
    print(f"Grid: {static_params['nx']}x{static_params['ny']}x{static_params['nz']}")
    print(f"Resolution: {static_params['dx']}m x {static_params['dy']}m x {static_params['dz']}m")
    print(f"Origin time: {static_params['origin_time']}")
    create_chemistry_driver(static_params)
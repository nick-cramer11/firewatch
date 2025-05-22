import cdsapi

# This script retrieves ERA5 data for a specific region and time period, processes it, and saves the output.
import os
import sys

os.environ["CDSAPI_RC"] = r"C:\Users\jginn\.cdsapirc"

if not os.path.exists(os.path.expanduser("~/.cdsapirc")):
    sys.exit("Error: CDS API key is not configured. Please ensure the '.cdsapirc' file exists in your home directory.")

c = cdsapi.Client()

c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            '2m_temperature', '10m_u_component_of_wind', '10m_v_component_of_wind',
            '2m_dewpoint_temperature', 'surface_pressure'
        ],
        'year': '2024',
        'month': ['04', '05', '10'],  # spring/fall months
        'day': ['01', '02', '03'],    # test subset
        'time': ['12:00'],
        'format': 'netcdf',
        'area': [41.5, -124.5, 39.5, -123.0],  # N, W, S, E — Humboldt Co. box
    },
    'era5_subset.nc')

import xarray as xr
import geopandas as gpd
import rioxarray

era5 = xr.open_dataset("era5_subset.nc")

# Convert temperature from Kelvin to Celsius
era5["t2m"] = era5["t2m"] - 273.15

# Clip using shapefile
shapefile_path = r"firewatch\data\hrc_ownership_polygon\hrc_ownership_final.shp"
gdf = gpd.read_file(shapefile_path)
gdf = gdf.to_crs("EPSG:4326")  # Convert to WGS 84
era5 = era5.rio.write_crs("EPSG:4326")  # Ensure the raster has the same CRS
era5_clip = era5.rio.clip(gdf.geometry.values, gdf.crs)

print("Shapefile Bounds:", gdf.total_bounds)
print("Raster Bounds:", era5.rio.bounds())

# Save output
output_path = r"firewatch\data\era5_humboldt.nc"
era5_clip.to_netcdf(output_path)
print(f"Saved clipped data to {output_path}")

# Subset the data for specific conditions
subset = era5_clip.where(
    (era5_clip["t2m"] > 10) & # temperature > 10 C
    (era5_clip["t2m"] < 25) & # and temperature < 25 C
    (((era5_clip["u10"]**2) + (era5_clip["v10"]**2))**0.5 < 10),  # and wind speed < 10 m/s
    drop=True # drop NaN values
)

# Save the subset
subset_path = r"firewatch\data\era5_humboldt_subset.nc"
subset.to_netcdf(subset_path)
print("ERA5 data subset saved to:", subset_path)
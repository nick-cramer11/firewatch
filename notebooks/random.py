import cdsapi

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
gdf = gpd.read_file("CNTYOUTL.SHP")
era5 = era5.rio.write_crs("EPSG:4326")
era5_clip = era5.rio.clip(gdf.geometry, gdf.crs)

# Save output
era5_clip.to_netcdf("era5_humboldt.nc")

subset = era5_clip.where(
    (era5_clip["t2m"] > 10) &
    (era5_clip["t2m"] < 25) &
    (era5_clip["u10"]**2 + era5_clip["v10"]**2)**0.5 < 10,  # wind speed < 10 m/s
    drop=True
)

import os
import sys
sys.path.append(r"C:\Users\jginn\OneDrive\Documents\OSU_24-25\spring25\GEOG562\project\python_code_project\firewatch\notebooks")
import era5_downloader
import xarray as xr
import netCDF4 as net
import geopandas as gpd
import importlib



importlib.reload(era5_downloader)

# retrieving ERA5 data for a specific year
year = 1990
output_dir = rf"firewatch/data/era5_data/{year}"
monthly_files = era5_downloader.retrieve_era5_data(year, output_dir = output_dir)
print(monthly_files)



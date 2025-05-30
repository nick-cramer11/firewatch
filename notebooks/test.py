import os
import sys
sys.path.append(r"C:\Users\jginn\OneDrive\Documents\OSU_24-25\spring25\GEOG562\project\python_code_project\firewatch\notebooks")
import data_import
import xarray as xr
import netCDF4 as net
import geopandas as gpd
import importlib



importlib.reload(data_import)

# retrieving ERA5 data for a specific year
year = 1990
output_dir = r"firewatch/data/era5_data/1990"
monthly_files = data_import.retrieve_era5_data(year, output_dir = output_dir)
print(monthly_files)


# # using clipping function to clip ERA5 data for a specific study area
# folder_path = r"firewatch/data/era5_data/2024"
# study_area_shp = "firewatch/data/hrc_ownership_polygon/hrc_ownership_final.shp"
# output_folder = r"firewatch/data/era5_data/clipped_2024"

# try:
#     data_import.clip_era5_data(folder_path, study_area_shp, output_folder)
#     print("Clipping completed successfully!")
# except Exception as e:
#     print(f"An error occurred during clipping: {e}")


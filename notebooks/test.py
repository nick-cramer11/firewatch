import sys
sys.path.append(r"C:\Users\jginn\OneDrive\Documents\OSU_24-25\spring25\GEOG562\project\python_code_project\firewatch\notebooks")
import xarray as xr
import geopandas as gpd
import data
print("It works!")

monthly_files = data.retrieve_era5_data(year = 2024)
print(monthly_files)
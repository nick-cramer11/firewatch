
# This script retrieves ERA5 data for a specific region and time period, processes it, and saves the output.
import os
import sys

os.environ["CDSAPI_RC"] = r"C:\Users\jginn\.cdsapirc"

if not os.path.exists(os.path.expanduser("~/.cdsapirc")):
    sys.exit("Error: CDS API key is not configured. Please ensure the '.cdsapirc' file exists in your home directory.")

import cdsapi

# Function to retrieve ERA5 data for a specific year and area
def retrieve_era5_data(year=2024, variables=None, area=None):
    if variables is None:
        variables = [
            "2m_dewpoint_temperature",
            "2m_temperature",
            "10m_u_component_of_wind",
            "10m_v_component_of_wind"
        ]
    
    if area is None:
        area = [41.5, -124.5, 39.5, -123]  # Default area from your example

    client = cdsapi.Client()
    monthly_data = {} # initialize a dictionary to store monthly data
    
    # Loop through each month
    for month in range(1, 13):  # 1 to 12 for January to December
        request = {
            "product_type": "reanalysis",
            "variable": variables,
            "year": str(year),
            "month": [f"{month:02d}"],  # Zero-padded month format -- e.g., "01" for January
            "day": [str(i).zfill(2) for i in range(1, 32)],  # Days formatted as zero-padded strings
            "daily_statistic": "daily_mean",
            "time_zone": "utc-08:00",
            "frequency": "6_hourly",
            "format": "netcdf",
            "area": area
        }
        
        output_file = f"era5_{month:02d}.nc" # e.g., "era5_01.nc" for January
        
        if os.path.exists(output_file): 
            print(f"File {output_file} already exists. Skipping download.")
            continue 

        client.retrieve("derived-era5-single-levels-daily-statistics", request, output_file)
        
        monthly_data[month] = output_file  # Store output file per month in the monthly data dictionary
        print(f"Downloaded data for {year}-{month:02d} to {output_file}")
    
    return monthly_data

# Example usage
monthly_files = retrieve_era5_data()
print(monthly_files)


import xarray as xr
import rioxarray as rxr
import geopandas as gpd

def clip_era5_to_study_area(monthly_files, study_area_shp):
    """
    Clips each month's ERA5 dataset to the specified study area.

    Parameters:
    - monthly_files: Dictionary mapping month (1-12) to ERA5 NetCDF file.
    - study_area_shp: Path to the study area shapefile.
    - output_prefix: Prefix for the output clipped files.

    Returns:
    - Dictionary mapping month to clipped NetCDF file name.
    """
    study_area = gpd.read_file(study_area_shp)  # Load study area geometry
    clipped_data = {} # initialize a dictionary to store clipped data

    for month, file_path in monthly_files.items():
        dataset = xr.open_dataset(file_path)  # Load ERA5 data
        dataset = dataset.rio.write_crs("EPSG:4326")  # Ensure proper CRS
        
        # Clip each variable using the study area
        clipped_dataset = dataset.rio.clip(study_area.geometry.values, study_area.crs) 
        
        # Save clipped data
        output_file = f"clipped_era5_{month:02d}.nc"
        clipped_dataset.to_netcdf(output_file)
        clipped_data[month] = output_file

    return clipped_data

# Example usage
study_area_path = "your_study_area.shp"
clipped_files = clip_era5_to_study_area(monthly_files, study_area_path)
print(clipped_files)

# This script retrieves ERA5 data for a specific region and time period, processes it, and saves the output.
### Note: This script is designed to be run in an environment where the CDS API is configured.
import os
import sys
import calendar
from datetime import datetime
import time

os.environ["CDSAPI_RC"] = r"C:\Users\jginn\.cdsapirc"

if not os.path.exists(os.path.expanduser("~/.cdsapirc")):
    sys.exit("Error: CDS API key is not configured. Please ensure the '.cdsapirc' file exists in your home directory.")

import cdsapi

def retrieve_era5_data(year=2024, variables=None, area=None, output_dir=".", retry_attempts=3):
    """
    Retrieve ERA5 data for a specific year and area.
    
    Args:
        year (int): Year to download data for
        variables (list): List of variables to download
        area (list): Bounding box [North, West, South, East]
        output_dir (str): Directory to save files
        retry_attempts (int): Number of retry attempts for failed downloads
    
    Returns:
        dict: Dictionary mapping months to their output files
    """
    if variables is None:
        variables = [
            "2m_dewpoint_temperature",
            "2m_temperature",
            "10m_u_component_of_wind",
            "10m_v_component_of_wind"
        ]
    
    if area is None:
        area = [41.5, -124.5, 39.5, -123]  # Default area

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    client = cdsapi.Client()
    monthly_data = {}
    failed_downloads = []
    
    print(f"Starting ERA5 download for year {year}")
    print(f"Variables: {', '.join(variables)}")
    print(f"Area: {area}")
    print("-" * 50)
    
    # Loop through each month
    for month in range(1, 13):
        # Get the correct number of days for this month/year
        days_in_month = calendar.monthrange(year, month)[1]
        days = [f"{i:02d}" for i in range(1, days_in_month + 1)]
        
        # Original working request parameters
        request = {
            "product_type": "reanalysis",
            "variable": variables,
            "year": str(year),
            "month": [f"{month:02d}"],  # Keep as list format
            "day": days,
            "daily_statistic": "daily_mean",
            "time_zone": "utc-08:00",  # Your original timezone
            "frequency": "6_hourly",   # Your original frequency
            "format": "netcdf",
            "area": area
        }
        
        output_file = os.path.join(output_dir, f"era5_{month:02d}_{year}.nc")
        
        if os.path.exists(output_file):
            print(f"✓ File {os.path.basename(output_file)} already exists. Skipping download.")
            monthly_data[month] = output_file
            continue
        
        # Attempt download with retry logic
        success = False
        for attempt in range(retry_attempts):
            try:
                print(f"📥 Downloading {calendar.month_name[month]} {year} (attempt {attempt + 1}/{retry_attempts})...")
                start_time = time.time()
                
                        # Use the correct dataset identifier for post-processed daily statistics
                client.retrieve("derived-era5-single-levels-daily-statistics", request, output_file)
                
                elapsed_time = time.time() - start_time
                file_size = os.path.getsize(output_file) / (1024 * 1024)  # Size in MB
                
                print(f"✓ Downloaded {os.path.basename(output_file)} ({file_size:.1f} MB) in {elapsed_time:.1f}s")
                monthly_data[month] = output_file
                success = True
                break
                
            except Exception as e:
                print(f"✗ Attempt {attempt + 1} failed: {str(e)}")
                if attempt < retry_attempts - 1:
                    wait_time = 2 ** attempt  # Exponential backoff
                    print(f"  Waiting {wait_time}s before retry...")
                    time.sleep(wait_time)
                else:
                    print(f"✗ Failed to download {calendar.month_name[month]} {year} after {retry_attempts} attempts")
                    failed_downloads.append((month, str(e)))
        
        if not success and os.path.exists(output_file):
            # Clean up partial download
            os.remove(output_file)
    
    # Summary
    print("-" * 50)
    print(f"Download Summary:")
    print(f"✓ Successfully downloaded: {len(monthly_data)}/12 months")
    
    if failed_downloads:
        print(f"✗ Failed downloads: {len(failed_downloads)}")
        for month, error in failed_downloads:
            print(f"  - {calendar.month_name[month]}: {error}")
    
    return monthly_data, failed_downloads if failed_downloads else monthly_data

def validate_parameters(year, variables, area):
    """Validate input parameters before making API requests."""
    current_year = datetime.now().year
    
    if not (1940 <= year <= current_year):
        raise ValueError(f"Year must be between 1940 and {current_year}")
    
    valid_variables = [
        "2m_dewpoint_temperature", "2m_temperature", "10m_u_component_of_wind",
        "10m_v_component_of_wind", "surface_pressure", "total_precipitation"
        # Add more as needed
    ]
    
    invalid_vars = [var for var in variables if var not in valid_variables]
    if invalid_vars:
        raise ValueError(f"Invalid variables: {invalid_vars}")
    
    if len(area) != 4:
        raise ValueError("Area must be [North, West, South, East]")
    
    north, west, south, east = area
    if not (-90 <= south < north <= 90):
        raise ValueError("Invalid latitude range")
    if not (-180 <= west < east <= 180):
        raise ValueError("Invalid longitude range")

# Example usage
if __name__ == "__main__":
    # Download with validation
    try:
        year = 2024
        variables = ["2m_temperature", "2m_dewpoint_temperature"]
        area = [41.5, -124.5, 39.5, -123]
        
        validate_parameters(year, variables, area)
        result = retrieve_era5_data(
            year=year,
            variables=variables,
            area=area,
            output_dir="./era5_data",
            retry_attempts=3
        )
        
        if isinstance(result, tuple):
            monthly_data, failed_downloads = result
            print(f"\nDownload completed with {len(failed_downloads)} failures")
        else:
            monthly_data = result
            print(f"\nAll downloads completed successfully!")
            
    except Exception as e:
        print(f"Error: {e}")



# This script retrieves ERA5 data for a specific region and time period, processes it, and saves the output.
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


# clipping function not working because of file permission issues with the NetCDF files
# going to skip this step for now and focus on retrieving the data
import os
import zipfile
import geopandas as gpd
import xarray as xr
from shapely.ops import unary_union
import tempfile
import shutil

def clip_era5_data(folder_path, study_area_shp, output_folder):
    """
    Uncompress ZIP files and clip all ERA5 NetCDF files in a folder using a shapefile.
    
    Parameters:
    -----------
    folder_path : str
        Path to folder containing NetCDF files (compressed or uncompressed)
    study_area_shp : str
        Path to shapefile defining the study area
    output_folder : str
        Path to output folder for clipped NetCDF files
    
    Returns:
    --------
    dict : Dictionary mapping original filenames to output file paths
    """
    
    # Load study area geometry
    study_area = gpd.read_file(study_area_shp)
    print(f"Loaded study area with CRS: {study_area.crs}")
    print(f"Study area shape: {study_area.shape}")
    
    # Ensure it's a Polygon/MultiPolygon
    study_area = study_area[study_area.geom_type.isin(["Polygon", "MultiPolygon"])]
    
    # Always reproject to WGS84 (EPSG:4326) since NetCDF files are in geographic coordinates
    print(f"Original study area bounds: {study_area.bounds.iloc[0].values}")
    if study_area.crs != 'EPSG:4326':
        print(f"Reprojecting study area from {study_area.crs} to EPSG:4326")
        study_area_wgs84 = study_area.to_crs('EPSG:4326')
        print(f"Reprojected bounds: {study_area_wgs84.bounds.iloc[0].values}")
    else:
        study_area_wgs84 = study_area
    
    # Merge geometries into a single object if needed
    geometry = unary_union(study_area_wgs84.geometry)
    print(f"Geometry type after union: {type(geometry)}")
    print(f"Final geometry bounds (WGS84): {geometry.bounds}")
    
    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)
    
    clipped_data = {}  # Store clipped files
    
    # Create temporary directory for extracted files
    temp_dir = tempfile.mkdtemp()
    
    try:
        # Get list of all files first
        all_files = os.listdir(folder_path)
        print(f"Found {len(all_files)} total files in folder: {all_files}")
        
        # Loop through all files in the folder
        for i, filename in enumerate(all_files):
            print(f"\n--- Processing file {i+1}/{len(all_files)}: {filename} ---")
            file_path = os.path.join(folder_path, filename)
            
            # Skip if it's a directory
            if not os.path.isfile(file_path):
                print(f"Skipping {filename} (not a file)")
                continue
            
            nc_files_to_process = []
            
            # Better file type detection
            if filename.endswith('.nc'):
                # Double-check it's actually a ZIP file
                if zipfile.is_zipfile(file_path):
                    print(f"Processing ZIP file: {filename}")
                    nc_files_to_process.extend(extract_zip_netcdf(file_path, temp_dir))
                else:
                    print(f"File {filename} has .zip extension but is not a valid ZIP file")
                    # Check if it might be a NetCDF file with wrong extension
                    try:
                        # Try to open as NetCDF
                        test_ds = xr.open_dataset(file_path, engine="netcdf4")
                        test_ds.close()
                        print(f"  -> Treating as NetCDF file instead")
                        nc_files_to_process.append(file_path)
                    except:
                        print(f"  -> Cannot open as NetCDF either, skipping")
                        continue
                
            elif filename.endswith(('.nc', '.netcdf', '.nc4')):
                print(f"Processing NetCDF file: {filename}")
                nc_files_to_process.append(file_path)
                
            else:
                print(f"Skipping {filename} (unsupported file type)")
                continue
            
            print(f"Found {len(nc_files_to_process)} NetCDF files to process from {filename}")
            
            # Process each NetCDF file
            for j, nc_file_path in enumerate(nc_files_to_process):
                print(f"  Processing NetCDF {j+1}/{len(nc_files_to_process)}: {os.path.basename(nc_file_path)}")
                try:
                    # Get original filename for output naming
                    original_name = os.path.basename(nc_file_path)
                    if nc_file_path.startswith(temp_dir):
                        # This came from a ZIP, use the ZIP name as prefix
                        zip_name = os.path.splitext(filename)[0]  # Remove .zip extension
                        original_name = f"{zip_name}_{original_name}"
                    
                    # Load NetCDF data
                    print(f"    Loading dataset: {original_name}")
                    dataset = xr.open_dataset(nc_file_path, engine="netcdf4")
                    print(f"    Dataset loaded successfully. Shape: {dataset.dims}")
                    print(f"    Dataset coordinates: {list(dataset.coords.keys())}")
                    
                    # Handle CRS - NetCDF data is typically in WGS84 but may not have CRS set
                    if dataset.rio.crs is None:
                        print(f"    Setting CRS to EPSG:4326 (WGS84)")
                        dataset = dataset.rio.write_crs("EPSG:4326")
                    else:
                        print(f"    Dataset CRS: {dataset.rio.crs}")
                    
                    # Get dataset bounds for debugging
                    try:
                        ds_bounds = dataset.rio.bounds()
                        print(f"    Dataset bounds: {ds_bounds}")
                    except Exception as e:
                        print(f"    Could not get dataset bounds: {e}")
                    
                    # Ensure coordinates are in the right order (some NetCDF files have lat/lon reversed)
                    coord_names = list(dataset.coords.keys())
                    lat_names = [name for name in coord_names if 'lat' in name.lower()]
                    lon_names = [name for name in coord_names if 'lon' in name.lower()]
                    
                    if lat_names and lon_names:
                        print(f"    Found latitude coord: {lat_names[0]}, longitude coord: {lon_names[0]}")
                        
                        # Check coordinate order and values
                        lat_vals = dataset[lat_names[0]].values
                        lon_vals = dataset[lon_names[0]].values
                        print(f"    Lat range: {lat_vals.min():.3f} to {lat_vals.max():.3f}")
                        print(f"    Lon range: {lon_vals.min():.3f} to {lon_vals.max():.3f}")
                    
                    # Convert geometry to the right format for clipping
                    if hasattr(geometry, '__iter__') and not hasattr(geometry, 'geom_type'):
                        geom_list = list(geometry)
                    else:
                        geom_list = [geometry]
                    
                    print(f"    Clipping with {len(geom_list)} geometry(ies)")
                    print(f"    Geometry bounds: {geometry.bounds}")
                    
                    # Clip entire dataset using the study area geometry
                    print(f"    Clipping dataset: {original_name}")
                    try:
                        # Use EPSG:4326 for clipping since both dataset and geometry are in WGS84
                        clipped_dataset = dataset.rio.clip(geom_list, "EPSG:4326", drop=True, invert=False)
                        print(f"    Dataset clipped successfully. New shape: {clipped_dataset.dims}")
                    except Exception as clip_error:
                        print(f"    Error during clipping: {clip_error}")
                        print(f"    Trying alternative clipping method...")
                        
                        # Alternative method: clip using bounding box first, then geometry
                        geom_bounds = geometry.bounds  # (minx, miny, maxx, maxy)
                        print(f"    Using bounding box: {geom_bounds}")
                        
                        try:
                            bbox_clipped = dataset.rio.clip_box(
                                minx=geom_bounds[0], miny=geom_bounds[1], 
                                maxx=geom_bounds[2], maxy=geom_bounds[3]
                            )
                            clipped_dataset = bbox_clipped.rio.clip(geom_list, "EPSG:4326", drop=True)
                            print(f"    Alternative clipping successful. New shape: {clipped_dataset.dims}")
                        except Exception as alt_error:
                            print(f"    Alternative clipping also failed: {alt_error}")
                            print(f"    Using bounding box clip only...")
                            clipped_dataset = dataset.rio.clip_box(
                                minx=geom_bounds[0], miny=geom_bounds[1], 
                                maxx=geom_bounds[2], maxy=geom_bounds[3]
                            )
                            print(f"    Bounding box clip successful. New shape: {clipped_dataset.dims}")
                    print(f"    Dataset clipped successfully")
                    
                    # Save clipped data
                    output_filename = f"clipped_{original_name}"
                    output_file = os.path.join(output_folder, output_filename)
                    print(f"    Saving to: {output_file}")
                    clipped_dataset.to_netcdf(output_file)
                    
                    # Store in results dictionary
                    clipped_data[original_name] = output_file
                    
                    print(f"    ✓ Successfully processed: {original_name}")
                    
                    # Close dataset to free memory
                    dataset.close()
                    clipped_dataset.close()
                    
                except Exception as e:
                    print(f"    ✗ Error processing {nc_file_path}: {e}")
                    print(f"    Error type: {type(e).__name__}")
                    import traceback
                    print(f"    Full traceback:\n{traceback.format_exc()}")
                    continue
                    
            # Clear the list for next iteration
            nc_files_to_process = []
    
    finally:
        # Clean up temporary directory
        print(f"\nCleaning up temporary directory: {temp_dir}")
        shutil.rmtree(temp_dir, ignore_errors=True)
    
    print(f"\n=== PROCESSING COMPLETE ===")
    print(f"Successfully processed {len(clipped_data)} files:")
    for original, clipped in clipped_data.items():
        print(f"  ✓ {original}")
    
    return clipped_data


def extract_zip_netcdf(zip_path, extract_dir):
    """
    Extract NetCDF files from a ZIP archive.
    
    Parameters:
    -----------
    zip_path : str
        Path to ZIP file
    extract_dir : str
        Directory to extract files to
    
    Returns:
    --------
    list : List of paths to extracted NetCDF files
    """
    nc_files = []
    
    try:
        # First check if it's actually a ZIP file
        if not zipfile.is_zipfile(zip_path):
            print(f"Warning: {os.path.basename(zip_path)} is not a valid ZIP file, treating as regular NetCDF")
            # If it has .nc extension, treat it as a regular NetCDF file
            if zip_path.endswith(('.nc', '.netcdf', '.nc4')):
                return [zip_path]
            else:
                return []
        
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Get list of files in the ZIP
            file_list = zip_ref.namelist()
            print(f"Files in ZIP {os.path.basename(zip_path)}: {file_list}")
            
            # Extract all files
            zip_ref.extractall(extract_dir)
            
            # Find NetCDF files in extracted content
            for file_name in file_list:
                if file_name.endswith(('.nc', '.netcdf', '.nc4')):
                    extracted_path = os.path.join(extract_dir, file_name)
                    if os.path.exists(extracted_path):
                        nc_files.append(extracted_path)
                        print(f"Found NetCDF file in ZIP: {file_name}")
            
            if not nc_files:
                print(f"Warning: No NetCDF files found in {os.path.basename(zip_path)}")
                
    except PermissionError as e:
        print(f"Permission error accessing {zip_path}: {e}")
        return []
    except Exception as e:
        print(f"Error extracting {zip_path}: {e}")
        return []
    
    return nc_files



# Example usage:
if __name__ == "__main__":
    folder_path = "path/to/your/netcdf/files"
    study_area_shp = "path/to/your/study_area.shp"
    output_folder = "path/to/output/folder"
    
    results = clip_era5_data(folder_path, study_area_shp, output_folder)
    print(f"Successfully processed {len(results)} files:")
    for original, clipped in results.items():
        print(f"  {original} -> {clipped}")
# function to calculate the number of days in a year with specific weather conditions
import xarray as xr
import numpy as np

def count_qualifying_days(nc_file, lat, lon, year, temp_range=(50, 90), humidity_range=(20, 80), wind_max=10):
    """
    Counts the number of days meeting temperature, humidity, and wind speed criteria for a given raster cell.
    
    Parameters:
        nc_file (str): Path to NetCDF file
        lat (float): Latitude of the cell
        lon (float): Longitude of the cell
        year (int): Year to analyze
        temp_range (tuple): (min_temp_F, max_temp_F)
        humidity_range (tuple): (min_humidity_%, max_humidity_%)
        wind_max (float): Maximum allowable wind speed (mph)

    Returns:
        int: Count of qualifying days
    """
    # Load dataset
    ds = xr.open_dataset(nc_file)
    
    # Select cell closest to lat/lon
    cell = ds.sel(latitude=lat, longitude=lon, method="nearest")

    # Convert temperature from K to F
    temp_F = (cell['t2m'] - 273.15) * 9/5 + 32
    dewpoint_F = (cell['d2m'] - 273.15) * 9/5 + 32
    
    # Compute relative humidity (simplified formula)
    rh = 100 * (np.exp((17.625 * dewpoint_F) / (243.04 + dewpoint_F)) /
                np.exp((17.625 * temp_F) / (243.04 + temp_F)))
    
    # Compute wind speed
    wind_speed = np.sqrt(cell['u10']**2 + cell['v10']**2) * 2.237  # Convert m/s to mph

    # Extract year and filter data
    time_index = cell['time'].dt.year == year
    qualifying_days = ((temp_F[time_index] >= temp_range[0]) & (temp_F[time_index] <= temp_range[1]) &
                       (rh[time_index] >= humidity_range[0]) & (rh[time_index] <= humidity_range[1]) &
                       (wind_speed[time_index] <= wind_max)).sum().item()

    return int(qualifying_days)
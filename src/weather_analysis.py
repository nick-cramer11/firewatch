
##################################
# function to calculate the number of days in a year with specific weather conditions at a specified location

import os
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime

def count_qualifying_days(data_folder, year, month, lat, lon, temp_range=(50, 90), humidity_range=(20, 80), wind_max=10):
    """Counts qualifying days using separate files for temperature, humidity, and wind."""

    # Construct file paths dynamically firewatch/data/era5_data/2024/2m_dewpoint_temperature_stream-oper_daily-mean_01_2024.nc
    dewpoint_file = os.path.join(data_folder, f"2m_dewpoint_temperature_stream-oper_daily-mean_{month:02d}_{year}.nc")
    temp_file = os.path.join(data_folder, f"2m_temperature_0_daily-mean_{month:02d}_{year}.nc")
    u_wind_file = os.path.join(data_folder, f"10m_u_component_of_wind_0_daily-mean_{month:02d}_{year}.nc")
    v_wind_file = os.path.join(data_folder, f"10m_v_component_of_wind_0_daily-mean_{month:02d}_{year}.nc")

    # Check if all files exist before proceeding
    missing_files = [f for f in [temp_file, dewpoint_file, u_wind_file, v_wind_file] if not os.path.exists(f)]
    if missing_files:
        print(f"Missing files for {month}/{year}: {missing_files}")
        return None

    # Open datasets
    temp_ds = xr.open_dataset(temp_file, engine="netcdf4")
    dewpoint_ds = xr.open_dataset(dewpoint_file, engine="netcdf4")
    u_wind_ds = xr.open_dataset(u_wind_file, engine="netcdf4")
    v_wind_ds = xr.open_dataset(v_wind_file, engine="netcdf4")

    # Select closest cell based on lat/lon
    temp_cell = temp_ds.sel(latitude=lat, longitude=lon, method="nearest")
    dewpoint_cell = dewpoint_ds.sel(latitude=lat, longitude=lon, method="nearest")
    u_wind_cell = u_wind_ds.sel(latitude=lat, longitude=lon, method="nearest")
    v_wind_cell = v_wind_ds.sel(latitude=lat, longitude=lon, method="nearest")

    # Convert temperature from K to F
    temp_F = (temp_cell['t2m'] - 273.15) * 9/5 + 32
    dewpoint_F = (dewpoint_cell['d2m'] - 273.15) * 9/5 + 32

    # Compute relative humidity
    rh = 100 * (np.exp((17.625 * dewpoint_F) / (243.04 + dewpoint_F)) /
                np.exp((17.625 * temp_F) / (243.04 + temp_F)))

    # Compute wind speed
    wind_speed = np.sqrt(u_wind_cell['u10']**2 + v_wind_cell['v10']**2) * 2.237  # Convert m/s to mph

    # Filter by month and year
    time_index = (temp_cell['valid_time'].dt.year == year) & (temp_cell['valid_time'].dt.month == month)
    
    # Calculate qualifying days
    qualifying_days = ((temp_F[time_index] >= temp_range[0]) & (temp_F[time_index] <= temp_range[1]) &
                       (rh[time_index] >= humidity_range[0]) & (rh[time_index] <= humidity_range[1]) &
                       (wind_speed[time_index] <= wind_max)).sum().item()

    return int(qualifying_days)

##################################
# function to count the number of qualifying days across all months in a year and create a raster-like output

def count_qualifying_days_raster(data_folder, year, month_range=(1, 12), temp_range=(50, 90), humidity_range=(20, 80), wind_max=10):
    """
    Counts the total number of qualifying days per raster cell across the entire year.

    Parameters:
        data_folder (str): Path to folder containing NetCDF files
        year (int): Year to analyze
        month_range (tuple): Range of months to analyze (e.g., (1, 12) for Jan-Dec)
        temp_range (tuple): (min_temp_F, max_temp_F)
        humidity_range (tuple): (min_humidity_%, max_humidity_%)
        wind_max (float): Maximum allowable wind speed (mph)

    Returns:
        xarray.DataArray: A raster-like array with the total qualifying day counts for each cell over the year.
    """
    # Initialize an empty raster to accumulate qualifying days across all months
    total_qualifying_days = None  

    for month in range(month_range[0], month_range[1] + 1):
        # Construct file paths dynamically
        dewpoint_file = os.path.join(data_folder, f"2m_dewpoint_temperature_stream-oper_daily-mean_{month:02d}_{year}.nc")
        temp_file = os.path.join(data_folder, f"2m_temperature_0_daily-mean_{month:02d}_{year}.nc")
        u_wind_file = os.path.join(data_folder, f"10m_u_component_of_wind_0_daily-mean_{month:02d}_{year}.nc")
        v_wind_file = os.path.join(data_folder, f"10m_v_component_of_wind_0_daily-mean_{month:02d}_{year}.nc")

        # Check if all files exist
        missing_files = [f for f in [temp_file, dewpoint_file, u_wind_file, v_wind_file] if not os.path.exists(f)]
        if missing_files:
            print(f"Skipping {month}/{year} due to missing files: {missing_files}")
            continue  # Skip this month, but keep processing other months

        # Load datasets
        temp_ds = xr.open_dataset(temp_file, engine="netcdf4")
        dewpoint_ds = xr.open_dataset(dewpoint_file, engine="netcdf4")
        u_wind_ds = xr.open_dataset(u_wind_file, engine="netcdf4")
        v_wind_ds = xr.open_dataset(v_wind_file, engine="netcdf4")

        # Convert 'valid_time' to datetime
        temp_ds['valid_time'] = temp_ds['valid_time'].astype("datetime64[ns]")

        # Convert temperature from K to F
        temp_F = (temp_ds['t2m'] - 273.15) * 9/5 + 32
        dewpoint_F = (dewpoint_ds['d2m'] - 273.15) * 9/5 + 32

        # Compute relative humidity
        rh = 100 * (np.exp((17.625 * dewpoint_F) / (243.04 + dewpoint_F)) /
                    np.exp((17.625 * temp_F) / (243.04 + temp_F)))

        # Compute wind speed
        wind_speed = np.sqrt(u_wind_ds['u10']**2 + v_wind_ds['v10']**2) * 2.237  # Convert m/s to mph

        # Filter by year and month
        time_index = (temp_ds['valid_time'].dt.year == year) & (temp_ds['valid_time'].dt.month == month)

        # Calculate qualifying days for the current month
        qualifying_days = ((temp_F[time_index] >= temp_range[0]) & (temp_F[time_index] <= temp_range[1]) &
                           (rh[time_index] >= humidity_range[0]) & (rh[time_index] <= humidity_range[1]) &
                           (wind_speed[time_index] <= wind_max)).sum(dim="valid_time")

        # Accumulate qualifying days across all months
        if total_qualifying_days is None:
            total_qualifying_days = qualifying_days  # Initialize the array with the first valid month
        else:
            total_qualifying_days += qualifying_days  # Add monthly counts to the total

    return total_qualifying_days

########################################

# function to create a raster where each cell displays the month with the most qualifying days

def find_best_month_raster(data_folder, year, month_range=(1, 12), temp_range=(50, 90), humidity_range=(20, 80), wind_max=10):
    """
    Finds the month with the highest number of qualifying days per raster cell.
    
    Parameters:
        data_folder (str): Path to folder containing NetCDF files
        year (int): Year to analyze
        month_range (tuple): Range of months to analyze (e.g., (1, 12) for Jan-Dec)
        temp_range (tuple): (min_temp_F, max_temp_F)
        humidity_range (tuple): (min_humidity_%, max_humidity_%)
        wind_max (float): Maximum allowable wind speed (mph)

    Returns:
        xarray.DataArray: Raster where each cell contains the **best month** (1–12) with the most qualifying days.
    """
    best_month_raster = None  # Stores the best month per cell
    max_days_raster = None    # Tracks the highest qualifying days per cell

    for month in range(month_range[0], month_range[1] + 1):
        # Construct file paths
        dewpoint_file = os.path.join(data_folder, f"2m_dewpoint_temperature_stream-oper_daily-mean_{month:02d}_{year}.nc")
        temp_file = os.path.join(data_folder, f"2m_temperature_0_daily-mean_{month:02d}_{year}.nc")
        u_wind_file = os.path.join(data_folder, f"10m_u_component_of_wind_0_daily-mean_{month:02d}_{year}.nc")
        v_wind_file = os.path.join(data_folder, f"10m_v_component_of_wind_0_daily-mean_{month:02d}_{year}.nc")

        # Skip missing months
        if not all(map(os.path.exists, [temp_file, dewpoint_file, u_wind_file, v_wind_file])):
            print(f"Skipping month {month}/{year}, missing files.")
            continue

        # Load datasets
        temp_ds = xr.open_dataset(temp_file)
        dewpoint_ds = xr.open_dataset(dewpoint_file)
        u_wind_ds = xr.open_dataset(u_wind_file)
        v_wind_ds = xr.open_dataset(v_wind_file)

        # Convert temperature from K to F
        temp_F = (temp_ds['t2m'] - 273.15) * 9/5 + 32
        dewpoint_F = (dewpoint_ds['d2m'] - 273.15) * 9/5 + 32

        # Compute relative humidity
        rh = 100 * (np.exp((17.625 * dewpoint_F) / (243.04 + dewpoint_F)) /
                    np.exp((17.625 * temp_F) / (243.04 + temp_F)))

        # Compute wind speed
        wind_speed = np.sqrt(u_wind_ds['u10']**2 + v_wind_ds['v10']**2) * 2.237  # Convert m/s to mph

        # Filter by month
        time_index = (temp_ds['valid_time'].dt.year == year) & (temp_ds['valid_time'].dt.month == month)
        qualifying_days = ((temp_F[time_index] >= temp_range[0]) & (temp_F[time_index] <= temp_range[1]) &
                           (rh[time_index] >= humidity_range[0]) & (rh[time_index] <= humidity_range[1]) &
                           (wind_speed[time_index] <= wind_max)).sum(dim="valid_time")

        # Initialize raster on first valid month
        if best_month_raster is None:
            best_month_raster = xr.full_like(qualifying_days, fill_value=month, dtype=int)
            max_days_raster = qualifying_days
        else:
            # Update only where this month's qualifying days exceed the previous max
            best_month_raster = xr.where(qualifying_days > max_days_raster, month, best_month_raster)
            max_days_raster = xr.where(qualifying_days > max_days_raster, qualifying_days, max_days_raster)

    return best_month_raster


# ========================================
# Script from all_months_all_years.py
# This script was moved for better organization.

# === Thresholds ===
TEMP_MIN_F = 50
TEMP_MAX_F = 90
TEMP_MIN = (TEMP_MIN_F - 32) * 5 / 9  # Convert to Celsius
TEMP_MAX = (TEMP_MAX_F - 32) * 5 / 9
RH_MIN = 20
RH_MAX = 80
WIND_MAX = 10  # m/s

# === Input folder with unzipped .nc files ===
unzipped_dir = "data/unzipped"

# === List all months by matching one of the files ===
months = sorted(list(set([
    "_".join(f.split("_")[-2:]).replace(".nc", "")  # e.g., "01_1990"
    for f in os.listdir(unzipped_dir)
    if f.startswith("2m_temperature_0_daily-mean")
])))

# === Initialize dataframe for HTML output ===
records = []

for month_year in months:
    try:
        print(f"🔍 Processing {month_year}...")
        temp_ds = xr.open_dataset(f"{unzipped_dir}/2m_temperature_0_daily-mean_{month_year}.nc")
        dew_ds = xr.open_dataset(f"{unzipped_dir}/2m_dewpoint_temperature_stream-oper_daily-mean_{month_year}.nc")
        u10_ds = xr.open_dataset(f"{unzipped_dir}/10m_u_component_of_wind_0_daily-mean_{month_year}.nc")
        v10_ds = xr.open_dataset(f"{unzipped_dir}/10m_v_component_of_wind_0_daily-mean_{month_year}.nc")

        t = temp_ds['t2m'] - 273.15  # Celsius
        td = dew_ds['d2m'] - 273.15
        u10 = u10_ds['u10']
        v10 = v10_ds['v10']

        rh = 100 * (np.exp((17.625 * td) / (243.04 + td)) / np.exp((17.625 * t) / (243.04 + t)))
        wind = np.sqrt(u10**2 + v10**2)

        for i in range(len(t['valid_time'])):
            t_day = t.isel(valid_time=i)
            rh_day = rh.isel(valid_time=i)
            wind_day = wind.isel(valid_time=i)

            t_mean = float(t_day.mean())
            rh_mean = float(rh_day.mean())
            wind_mean = float(wind_day.mean())

            date = pd.to_datetime(t['valid_time'][i].values).strftime("%Y-%m-%d")
            is_eligible = (
                TEMP_MIN <= t_mean <= TEMP_MAX and
                RH_MIN <= rh_mean <= RH_MAX and
                wind_mean <= WIND_MAX
            )
            symbol = "Yes" if is_eligible else "No"

            records.append({
                "Date": date,
                "Month_Year": month_year,
                "Temp_C": round(t_mean, 2),
                "RH_%": round(rh_mean, 2),
                "Wind_m/s": round(wind_mean, 2),
                "Eligible": symbol
            })

    except Exception as e:
        print(f"⚠️ Failed to process {month_year}: {e}")

# === Save as HTML table ===
df = pd.DataFrame(records)
df.to_html("burn_day_report.html", index=False, escape=False)
print("✅ HTML report saved as 'burn_day_report.html2'")

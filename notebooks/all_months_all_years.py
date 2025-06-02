import os
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime

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
            symbol = "✅" if is_eligible else "❌"

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

import xarray as xr
import numpy as np
import os

data_root = "data/unzipped"


def count_qualifying_days_raster(data_folder, year):
    try:
        # Thresholds
        TEMP_MIN = (50 - 32) * 5/9
        TEMP_MAX = (90 - 32) * 5/9
        RH_MIN = 20
        RH_MAX = 80
        WIND_MAX = 10  # m/s

        qualifying_sum = None

        # Loop through all 12 months
        for month in range(1, 13):
            suffix = f"{month:02d}_{year}"
            try:
                # Load datasets for this month
                t2m = xr.open_dataset(f"{data_folder}/2m_temperature_0_daily-mean_{suffix}.nc")['t2m'] - 273.15
                d2m = xr.open_dataset(f"{data_folder}/2m_dewpoint_temperature_stream-oper_daily-mean_{suffix}.nc")['d2m'] - 273.15
                u10 = xr.open_dataset(f"{data_folder}/10m_u_component_of_wind_0_daily-mean_{suffix}.nc")['u10']
                v10 = xr.open_dataset(f"{data_folder}/10m_v_component_of_wind_0_daily-mean_{suffix}.nc")['v10']
                
                # Calculate RH and wind speed
                rh = 100 * (np.exp((17.625 * d2m) / (243.04 + d2m)) / np.exp((17.625 * t2m) / (243.04 + t2m)))
                wind = np.sqrt(u10**2 + v10**2)

                # Burn day criteria
                mask = (
                    (t2m >= TEMP_MIN) & (t2m <= TEMP_MAX) &
                    (rh >= RH_MIN) & (rh <= RH_MAX) &
                    (wind <= WIND_MAX)
                )

                # Count qualifying days for this month
                count = mask.sum(dim="valid_time")

                # Accumulate across months
                if qualifying_sum is None:
                    qualifying_sum = count
                else:
                    qualifying_sum += count

            except FileNotFoundError:
                print(f"⚠️ Missing data for {suffix}, skipping...")
            except Exception as e:
                print(f"⚠️ Error processing {suffix}: {e}")

        return qualifying_sum

    except Exception as e:
        print(f"⚠️ Failed to process year {year}: {e}")
        return None

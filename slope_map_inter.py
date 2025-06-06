import os
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import weather as dpy

# === Setup ===
data_root = "data/unzipped"  # path to the folder with .nc files
intervals = [
    [1990, 1995, 2000],
    [2000, 2005, 2010],
    [2010, 2015, 2020]
]

# === Create slope maps for each interval ===
for period in intervals:
    raster_stack = []
    valid_years = []

    for year in period:
        print(f"Processing year: {year}")
        raster = dpy.count_qualifying_days_raster(data_root, year)
        if raster is not None:
            raster_stack.append(raster)
            valid_years.append(year)
        else:
            print(f"⚠️ Skipping year {year} due to missing data.")

    if len(raster_stack) < 2:
        print(f"⚠️ Not enough data to compute slope for {period}, skipping.")
        continue

    # Combine and assign coordinates
    combined = xr.concat(raster_stack, dim="year")
    combined = combined.assign_coords(year=("year", valid_years))

    # Compute slope
    slope = xr.apply_ufunc(
        lambda y: np.polyfit(valid_years, y, deg=1)[0],
        combined,
        input_core_dims=[["year"]],
        vectorize=True,
        output_dtypes=[float]
    )

    # Plot and save
    plt.figure(figsize=(10, 6))
    slope.plot(cmap="RdBu", vmin=-1, vmax=1)
    plt.title(f"Slope of Burn Days ({valid_years[0]}–{valid_years[-1]})")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.savefig(f"slope_map_{valid_years[0]}_{valid_years[-1]}.png", dpi=300)
    plt.close()

    print(f"✅ Saved: slope_map_{valid_years[0]}_{valid_years[-1]}.png")

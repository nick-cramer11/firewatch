import os
import xarray as xr
import numpy as np
import weather as dpy
import matplotlib.pyplot as plt


years = [1990, 1995, 2000, 2005, 2010, 2015, 2020, 2024]
data_root = "data/unzipped"

raster_stack = []
valid_years = []

for year in years:
    print(f"Processing year: {year}")
    raster = dpy.count_qualifying_days_raster(data_root, year)
    if raster is not None:
        raster_stack.append(raster)
        valid_years.append(year)
    else:
        print(f"⚠️ Skipping year {year} due to missing or incomplete data.")

if not raster_stack:
    print("No valid data found for any of the specified years. Exiting.")
    exit()

combined = xr.concat(raster_stack, dim="year")
combined = combined.assign_coords(year=("year", valid_years))



# Apply regression using polyfit
slope = xr.apply_ufunc(
    lambda y: np.polyfit(valid_years, y, deg=1)[0],  # slope only
    combined,
    input_core_dims=[["year"]],
    vectorize=True,
    dask="parallelized",
    output_dtypes=[float]
)


plt.figure(figsize=(10, 6))
slope.plot(cmap="RdBu", vmin=-1, vmax=1)  # xarray handles colorbar internally
plt.title("Slope of Eligible Burn Days for 1990–2024")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.savefig("slope_map.png", dpi=300)
plt.show()

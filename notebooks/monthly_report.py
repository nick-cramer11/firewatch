import os
import sys
sys.path.append(r"C:\Users\jginn\OneDrive\Documents\OSU_24-25\spring25\GEOG562\project\python_code_project\firewatch\notebooks")
# Import the custom module for counting qualifying days
import days_per_year as dpy
import unzip_and_rename as unzip
import importlib
importlib.reload(dpy)
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import pandas as pd

# make sure the working directory is set to the project root
os.chdir(r"C:\Users\jginn\OneDrive\Documents\OSU_24-25\spring25\GEOG562\project\python_code_project")

##############################
importlib.reload(dpy)
# unzip the data files if not already unzipped
year = 2015
unzip.unzip_and_rename_all(rf"firewatch/data/era5_data/{year}", rf"firewatch/data/era5_data/{year}")

#################################
### Calculating monthly stats

importlib.reload(dpy)
# Count qualifying days for each month in {year} at a specific location
data_folder = rf"firewatch/data/era5_data/{year}" # Define the directory containing your NetCDF files
lat = 40.606546
lon = -124.038802
monthly_results = {}

for month in range(1, 13):
    qualifying_days = dpy.count_qualifying_days(data_folder, year, month, lat, lon)
    monthly_results[month] = qualifying_days if qualifying_days is not None else "Data Missing"

# Print the results
print(monthly_results)

# Calculate the total qualifying days for the year
total_qualifying_days = sum(monthly_results.values())

for month in range(1, 13):
    # Calculate the percentage for one month (month 1)
    if total_qualifying_days > 0:
        month_percentage = (monthly_results[month] / total_qualifying_days) * 100
    else:
        month_percentage = 0  # Avoid division by zero if no data exists

    print(f"Month {month} contributed {month_percentage:.2f}% of the qualifying days in {year}.")

# calculate the percentage of qualifying days for each month
for month in range(1, 13):
    # Calculate the percentage for one month (month 1)
    if total_qualifying_days > 0:
        month_total = monthly_results[month]
    else:
        month_total = 0  # Avoid division by zero if no data exists

    print(f"Month {month} contributed {month_total} of the qualifying days in {year}.")


##########################################

### Plotting the results for the monthly qualifying days percentages

importlib.reload(dpy)
# Extract percentages for plotting
months = list(range(1, 13))
percentages = [(monthly_results[m] / total_qualifying_days) * 100 if total_qualifying_days > 0 else 0 for m in months]

# Create the plot
plt.figure(figsize=(10, 5))
plt.plot(months, percentages, marker="o", linestyle="-", color="royalblue")
# label the graph
plt.xlabel("Month")
plt.ylabel("Percentage of Yearly Qualifying Days")
plt.title(f"Percentage of Qualifying Days per Month in {year}")
plt.xticks(months)  # Ensure months are labeled correctly
plt.ylim(0, max(percentages) + 5)  # Set y-axis to fit values nicely

# Show the plot
plt.show()

#################################################

### Saving the results as an HTML table and generating a plot

importlib.reload(dpy)
# Extract percentages and totals for each month
percentages = [(monthly_results[m] / total_qualifying_days) * 100 if total_qualifying_days > 0 else 0 for m in months]
totals = [monthly_results[m] if total_qualifying_days > 0 else 0 for m in months]

# Create data frames for easy HTML conversion
data = {"Month": list(range(1, 13)), "Percentage": percentages}
df = pd.DataFrame(data)
df["Percentage"] = df["Percentage"].round(2)  # Round percentages to 2 decimal places

data2 = {"Month": list(range(1, 13)), "Qualifying Days": monthly_results.values()}
df2 = pd.DataFrame(data2)

# Save results as an HTML table
html_output = df.to_html(index=False)
html_output2 = df2.to_html(index=False)

# Generate a plot and save it as an image
plt.figure(figsize=(10, 5))
plt.plot(months, percentages, marker="o", linestyle="-", color="royalblue")
plt.xlabel("Month")
plt.ylabel("Percentage of Yearly Qualifying Days")
plt.title(f"Percentage of Qualifying Days per Month in {year}")
plt.xticks(df["Month"])
plt.ylim(0, max(df["Percentage"]) + 5)
plt.savefig(f"qualifying_days_plot_{year}.png")  # Save as image

# Combine results into an HTML file
html_report = f"""
<html>
<head>
    <title>Monthly Report {year}</title>
</head>
<body>
    <h1>Monthly Report for {year}</h1> <!-- Visible Title -->
    <h2>Qualifying Days per Month</h2>
    {html_output2}
    <br>
    <h2>Percentage of Yearly Qualifying Days per Month</h2>
    {html_output}
    <br><img src="qualifying_days_plot_{year}.png" width="1000">
</body>
</html>
"""

# Save HTML file
with open(f"qualifying_days_report_{year}.html", "w") as f:
    f.write(html_report)

print(f"HTML report saved: qualifying_days_report_{year}.html")

#######################################

### Raster of monthly qualifying days
importlib.reload(dpy)

months = (3, 6)  # Example: March to June
qualifying_days_raster = dpy.count_qualifying_days_raster(data_folder, year)

plt.figure(figsize=(8, 6))
qualifying_days_raster.plot(cmap="viridis")  # Change colormap if needed
plt.title(f"Qualifying Days Raster for {year}")
plt.xlabel("Longitude")
plt.ylabel("Latitude")
plt.colorbar(label="Number of Qualifying Days")  
plt.show()

#############################
### Raster of months with highest number of qualifying days
importlib.reload(dpy)

best_month_raster = dpy.find_best_month_raster(data_folder, year)
print(best_month_raster)

# Define custom colors for months
month_colors = {
    1: "red", 2: "orange", 3: "yellow", 4: "lightgreen", 5: "green",
    6: "cyan", 7: "blue", 8: "darkblue", 9: "purple", 10: "magenta",
    11: "pink", 12: "brown"
}
cmap = mcolors.ListedColormap([month_colors[m] for m in range(1, 13)])

# Create a figure and plot the raster
plt.figure(figsize=(8, 6))
best_month_raster.plot(cmap=cmap, add_colorbar=False)
plt.title(f"Month with Highest Qualifying Days in {year} in Humboldt County")
plt.xlabel("Longitude")
plt.ylabel("Latitude")

# Create manual legend
legend_patches = [mpatches.Patch(color=color, label=f"Month {month}") for month, color in month_colors.items()]
plt.legend(handles=legend_patches, bbox_to_anchor=(1.02, .5), loc='center left', title="Months")
plt.subplots_adjust(right=0.75)
plt.show()
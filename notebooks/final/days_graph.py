import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# Load and filter data
df = pd.read_html("burn_day_report.html")[0]
df = df[df["Eligible"] == "Yes"]
df["Year"] = pd.to_datetime(df["Date"]).dt.year
yearly_counts = df.groupby("Year").size().reset_index(name="Eligible_Burn_Days")

# Linear regression
slope, intercept, r_value, p_value, std_err = linregress(
    yearly_counts["Year"], yearly_counts["Eligible_Burn_Days"]
)
yearly_counts["Trend"] = intercept + slope * yearly_counts["Year"]

# Plot
plt.figure(figsize=(12, 6))

# Uniform point color (blue)
plt.scatter(
    yearly_counts["Year"],
    yearly_counts["Eligible_Burn_Days"],
    color="steelblue",
    s=100,
    edgecolor='black',
    zorder=3
)

# Connect points with dashed line
plt.plot(
    yearly_counts["Year"],
    yearly_counts["Eligible_Burn_Days"],
    color='grey',
    alpha=0.5,
    linestyle='--',
    zorder=2
)

# Trend line
plt.plot(
    yearly_counts["Year"],
    yearly_counts["Trend"],
    color="blue",
    linestyle="solid",
    label=f"Slope = {slope:.2f} days/year",
    zorder=1
)

# Add value labels above each point
for _, row in yearly_counts.iterrows():
    plt.text(
        row["Year"],
        row["Eligible_Burn_Days"] + 1,
        str(row["Eligible_Burn_Days"]),
        ha='center',
        fontsize=8
    )

# Labels and aesthetics
plt.title("Annual Eligible Burn Days Over Time", fontsize=16)
plt.xlabel("Year")
plt.ylabel("Number of Eligible Burn Days")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()

# Save as PNG
plt.savefig("burn_days_trend.png", dpi=300)

# Show the plot
plt.show()

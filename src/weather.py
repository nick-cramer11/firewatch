import xarray as xr
import json

def load_params(config_path):
    with open(config_path) as f:
        return json.load(f)

def filter_weather(ds, params):
    ds["t2m"] = ds["t2m"] - 273.15
    wind_speed = (ds["u10"]**2 + ds["v10"]**2)**0.5

    return ds.where(
        (ds["t2m"] > params["temperature_min"]) &
        (ds["t2m"] < params["temperature_max"]) &
        (wind_speed < params["wind_speed_max"]),
        drop=True
    )

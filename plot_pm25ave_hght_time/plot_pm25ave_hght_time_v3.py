import xarray as xr
import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
import os

# Define the file pattern
file_pattern = "Data/diagb_tracer/bkg_fv_tracer.202009??.*.nc"
file_list = sorted(glob.glob(file_pattern))

# Lists to hold results
pm25_profiles = []
times = []

# Parse datetime and compute spatial averages
for fname in file_list:
    try:
        ds = xr.open_dataset(fname)
        pm25 = ds['pm25_tot'].squeeze()
        pm25_avg_xy = pm25.mean(dim=['yaxis_1', 'xaxis_1'])
        pm25_profiles.append(pm25_avg_xy.values)

        # Extract filename only (not full path) for safer split
        base = os.path.basename(fname)
        parts = base.split('.')
        date_str = parts[1]  # e.g., '20200921'
        hour_str = parts[2]  # e.g., '06'
        time = pd.to_datetime(date_str + hour_str, format='%Y%m%d%H')
        times.append(time)
    except Exception as e:
        print(f"Error processing {fname}: {e}")

# Convert results into a DataArray
pm25_profiles = np.stack(pm25_profiles)  # shape: [time, z]
height_dim = ds['zaxis_1'].values

pm25_da = xr.DataArray(
    pm25_profiles,
    coords={'time': times, 'z': height_dim},
    dims=['time', 'z'],
    name='pm25_mean_profile'
)

# === Save the time–height profile to a NetCDF file ===
output_ds = pm25_da.to_dataset(name='pm25_mean_profile')
output_ds.to_netcdf("pm25_time_height.nc")
print("Saved time–height PM2.5 profile to pm25_time_height.nc")

# Plotting: time-height colormap
plt.figure(figsize=(10, 6))
pm25_da.T.plot(
    x='time',
    y='z',
    cmap='viridis',
    cbar_kwargs={'label': 'Mean PM₂.₅ (µg/m³)'},
    yincrease=False
)
plt.title("Time–Height PM₂.₅ Mean Profile (Spatially Averaged)")
plt.xlabel("Time (UTC)")
plt.ylabel("Model Level (z)")
plt.tight_layout()
plt.savefig("pm25_time_height.png", dpi=300)
plt.show()


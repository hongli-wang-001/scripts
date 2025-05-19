import os
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import matplotlib.dates as mdates
import csv

# Define experiment directories
experiment_dirs = [
    "AirNow",
    "vrf_3hcyc_ctr",
    "vrf_3hcyc_pm_h100v12",
    "vrf_3hcyc_pm_h100v5"
]
init_utc = "12"
#init_utc = "??"

var_obs_name = 'pm25'
var_fcst_name = 'pm25_diff'
exp_labels = []

# Data storage
exp_data = {exp: {'time': [], var_obs_name: [], var_fcst_name: [],'forecast_pm25': []} for exp in experiment_dirs}

# File search and data collection
for exp_index, exp_dir in enumerate(experiment_dirs):
    exp_name = os.path.basename(exp_dir)
    exp_labels.append(exp_name)

    file_pattern = os.path.join(exp_dir, f"2020*{init_utc}", "pm25.2020*_f???.nc")
    #file_pattern = os.path.join(exp_dir, f"2020*{init_utc}", "pm25.2020*_f??[1-6].nc")
    file_list = sorted(glob.glob(file_pattern))

    for file_path in file_list:
        try:
            # Extract datetime from file name
            basename = os.path.basename(file_path)
            datetime_str = basename.split('.')[1].split('_')[0]  # '2020092208'
            file_time = datetime.strptime(datetime_str, "%Y%m%d%H")

            # Load NetCDF file
            ds = xr.open_dataset(file_path)

            # Check if the required variables exist in the dataset
            if var_obs_name in ds and var_fcst_name in ds:
                # Calculate means, skipping NaN values
                pm25_mean = float(ds[var_obs_name].mean(skipna=True).values)
                pm25_diff_mean = float(ds[var_fcst_name].mean(skipna=True).values)

                # Forecast = obs + diff for experiments other than 'AirNow'
                if exp_index == 0:  # Assuming AirNow is the observed data
                    forecast_pm25_mean = pm25_mean #float(ds[var_obs_name].mean(skipna=True).values)
                else:
                    forecast_pm25_mean = pm25_mean + pm25_diff_mean # float(ds[var_fcst_name].mean(skipna=True).values)
                    #forecast_pm25_mean = forecast_pm25_mean + pm25_mean  # For others, we calculate forecast

                # Store data
                exp_data[exp_dir]['time'].append(file_time)
                exp_data[exp_dir][var_obs_name].append(pm25_mean)
                exp_data[exp_dir][var_fcst_name].append(pm25_diff_mean)
                exp_data[exp_dir]['forecast_pm25'].append(forecast_pm25_mean)
            else:
                # If the required variables are missing, append NaN values
                exp_data[exp_dir]['time'].append(file_time)
                exp_data[exp_dir][var_obs_name].append(np.nan)
                exp_data[exp_dir][var_fcst_name].append(np.nan)
                exp_data[exp_dir]['forecast_pm25'].append(np.nan)

            ds.close()

        except Exception as e:
            print(f"Error processing file {file_path}: {e}")


# --- Write forecast PM2.5 data to CSV ---
csv_filename = f"pm25_meanoverdom_obs_aqm_multiday_i{init_utc}.csv"

# Use the first experiment's time array as the reference
reference_times = exp_data[experiment_dirs[0]]['time']

# Open CSV file for writing
with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    
    # Write header row
    header = ['Time'] + exp_labels
    writer.writerow(header)

    # Write data row-by-row (one row per time)
    for i, t in enumerate(reference_times):
        row = [t.strftime('%Y-%m-%d %H:%M')]
        for exp_dir in experiment_dirs:
            value = exp_data[exp_dir]['forecast_pm25'][i]
            row.append(f"{value:.3f}" if not np.isnan(value) else "")
        writer.writerow(row)

print(f"Forecast PM2.5 time series written to: {csv_filename}")

# Plotting
fig, axes = plt.subplots(1, 1, figsize=(8, 4), sharex=True)

colors = ['black', 'blue', 'green', 'red']  # One color per experiment
#cmap = plt.get_cmap('tab10')  # You can use 'Set1', 'Dark2', etc.
#colors = [cmap(i) for i in range(len(experiment_dirs))]

# Plot Forecast PM2.5 = Obs + Diff
for idx, exp_dir in enumerate(experiment_dirs):
    axes.plot(exp_data[exp_dir]['time'], exp_data[exp_dir]['forecast_pm25'], label=exp_dir, color=colors[idx])

axes.set_title("Time Series of 1PM2.5 Mean")
axes.set_ylabel('PM2.5 (µg/m³)')
axes.set_xlabel("Time")
axes.legend()

# Improve x-axis: set range and add minor ticks
all_times = []
for exp_dir in experiment_dirs:
    all_times.extend(exp_data[exp_dir]['time'])
if all_times:
    time_min = min(all_times)
    time_max = max(all_times)
    axes.set_xlim([time_min, time_max])
    #axes.xaxis.set_major_locator(mdates.AutoDateLocator())
    axes.xaxis.set_major_locator(mdates.HourLocator(interval=24))
    axes.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d:%H'))
    #axes.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d\n%H:%M'))
    #axes.xaxis.set_minor_locator(mdates.HourLocator(interval=1))
    axes.tick_params(axis='x', which='major', rotation=60, length=8)
    #axes.tick_params(axis='x', which='minor', length=4, labelsize=8)

plt.tight_layout()
plt.savefig(f"pm25_meanoverdom_obs_aqm_multiday_i{init_utc}.png", dpi=300)
plt.show()

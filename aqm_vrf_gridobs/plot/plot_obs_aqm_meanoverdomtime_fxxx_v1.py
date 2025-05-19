import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob
import os
from collections import defaultdict
import csv

# List of experiment directories
experiment_dirs = [
    "AirNow",
    "vrf_3hcyc_ctr",
    "vrf_3hcyc_pm_h100v12",
    "vrf_3hcyc_pm_h100v5"
]
init_utc = "12"
forecast_hours = ['f001', 'f002', 'f003', 'f004', 'f005', 'f006', 'f007', 'f008',
                  'f009', 'f010', 'f011', 'f012', 'f013', 'f014', 'f015', 'f016',
                  'f017', 'f018', 'f019', 'f020', 'f021', 'f023', 'f024']
data_by_fhr = defaultdict(list)  # maps each fhr to list of data arrays (one per experiment)
exp_labels = []

var_obs_name = 'pm25_obs__mean'
# --- Data collection ---
#for exp_dir in experiment_dirs:
for exp_index, exp_dir in enumerate(experiment_dirs):
    exp_name = os.path.basename(exp_dir)
    exp_labels.append(exp_name)

    for fhr in forecast_hours:
        file_pattern = os.path.join(exp_dir, f"i{init_utc}", f"pm25_2d_stats_{fhr}.nc")
        files = sorted(glob.glob(file_pattern))
        all_data = []

        for file in files:
            with nc.Dataset(file, 'r') as ds:
                #if 'pm25_obs__mean' not in ds.variables:
                #    continue
                #data = ds.variables['pm25_obs__mean'][:]
                # Choose variable name based on experiment index
                var_name = 'pm25_obs__mean' if exp_index == 0 else 'pm25_diff_mean'
                if var_name not in ds.variables:
                    continue
                data = ds.variables[var_name][:]
                data_flat = data.flatten()
                data_clean = data_flat[~np.isnan(data_flat)]
                print(f"{os.path.basename(file)}: data_clean size = {data_clean.size}")
                if exp_index != 0:
                   data_obs = ds.variables[var_obs_name][:]
                   data_obs_flat = data_obs.flatten()
                   data_obs_clean = data_obs_flat[~np.isnan(data_obs_flat)]
                   print(f"{os.path.basename(file)}: data_obs_clean size = {data_obs_clean.size}")
                   data_clean = data_clean + data_obs_clean
 
                if data_clean.size > 0:
                    all_data.append(data_clean)

        # Combine and store
        if all_data:
            combined = np.concatenate(all_data)
            data_by_fhr[fhr].append(combined)
        else:
            data_by_fhr[fhr].append(np.array([]))  # placeholder for empty data

# --- Prepare data for line plot ---
mean_values = [[] for _ in range(len(experiment_dirs))]  # one list per experiment

for fhr in forecast_hours:
    for exp_index in range(len(experiment_dirs)):
        data = data_by_fhr[fhr][exp_index]
        if data.size > 0:
            mean_val = np.mean(data)
        else:
            mean_val = np.nan  # handle missing data
        mean_values[exp_index].append(mean_val)

# --- Save to CSV ---
csv_filename = f"pm25_meanoverdomtime_obs_aqm_i{init_utc}.csv"
with open(csv_filename, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    
    # Write header
    header = ["forecast_hour"] + exp_labels
    writer.writerow(header)
    
    # Write rows: one per forecast hour
    for i, fhr in enumerate(forecast_hours):
        row = [fhr] + [f"{mean_values[j][i]:.3f}" for j in range(len(experiment_dirs))]
        writer.writerow(row)

print(f"Mean values written to {csv_filename}")
# --- Plotting ---
colors = ['black', '#0072B2', '#E69F00', '#009E73']  # blue, orange, green
plt.figure(figsize=(10, 6))

for i, means in enumerate(mean_values):
    plt.plot(forecast_hours, means, label=exp_labels[i], color=colors[i], marker='o')

plt.xlabel("Forecast Hour")
plt.ylabel("Mean PM₂.₅ (µg/m³)")
plt.title("Mean PM₂.₅ by Forecast Hour")
plt.grid(True)
plt.legend(title='Experiment')
plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"pm25_meanoverdomtime_obs_aqm_i{init_utc}.png", dpi=300)
plt.show()


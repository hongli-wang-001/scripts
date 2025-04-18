import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob

# List of NetCDF files
file_pattern = './3dvar_diffusion_pm25_*.nc'
files = sorted(glob.glob(file_pattern))

# Initialize lists to store statistics
timestamps = []
obs_means, obs_medians, obs_stds = [], [], []
hofx0_means, hofx0_medians, hofx0_stds = [], [], []
hofx1_means, hofx1_medians, hofx1_stds = [], [], []

# Process each file
for file in files:
    with Dataset(file, 'r') as nc:
        # Extract time from the filename
        timestamp = file.split('_')[-1].split('.')[0]  # e.g., 2020090118
        timestamps.append(timestamp)

        # Extract the variables (update these names if needed)
        obs_var_name = 'ObsValue/particulatematter2p5Surface'  # Update this
        hofx0_var_name = 'hofx0/particulatematter2p5Surface'  # Update this
        hofx1_var_name = 'hofx1/particulatematter2p5Surface'  # Update this

        # Initialize lists for the current file
        obs_values = []
        hofx0_values = []
        hofx1_values = []

        # Check for variable existence and extract values
        if obs_var_name in nc.variables:
            obs_values = nc.variables[obs_var_name][:]
        if hofx0_var_name in nc.variables:
            hofx0_values = nc.variables[hofx0_var_name][:]
        if hofx1_var_name in nc.variables:
            hofx1_values = nc.variables[hofx1_var_name][:]

        # Calculate mean, median, and std for each variable
        if obs_values.size > 0:
            obs_means.append(np.mean(obs_values))
            obs_medians.append(np.median(obs_values))
            obs_stds.append(np.std(obs_values))
        else:
            obs_means.append(np.nan)
            obs_medians.append(np.nan)
            obs_stds.append(np.nan)

        if hofx0_values.size > 0:
            hofx0_means.append(np.mean(hofx0_values))
            hofx0_medians.append(np.median(hofx0_values))
            hofx0_stds.append(np.std(hofx0_values))
        else:
            hofx0_means.append(np.nan)
            hofx0_medians.append(np.nan)
            hofx0_stds.append(np.nan)

        if hofx1_values.size > 0:
            hofx1_means.append(np.mean(hofx1_values))
            hofx1_medians.append(np.median(hofx1_values))
            hofx1_stds.append(np.std(hofx1_values))
        else:
            hofx1_means.append(np.nan)
            hofx1_medians.append(np.nan)
            hofx1_stds.append(np.nan)

# Convert timestamps to a range suitable for plotting
x_ticks = range(len(timestamps))

# Create plots for means, medians, and standard deviations
plt.figure(figsize=(15, 10))

# Mean values
plt.subplot(3, 1, 1)
plt.plot(x_ticks, obs_means, label='ObsValue Mean', color='blue', marker='o')
plt.plot(x_ticks, hofx0_means, label='hofx0 Mean', color='green', marker='o')
plt.plot(x_ticks, hofx1_means, label='hofx1 Mean', color='purple', marker='o')
plt.title('Mean Values of Particulate Matter')
plt.xticks(x_ticks, timestamps, rotation=45)
plt.ylabel('Mean PM2.5')
plt.legend()

# Median values
plt.subplot(3, 1, 2)
plt.plot(x_ticks, obs_medians, label='ObsValue Median', color='blue', marker='o')
plt.plot(x_ticks, hofx0_medians, label='hofx0 Median', color='green', marker='o')
plt.plot(x_ticks, hofx1_medians, label='hofx1 Median', color='purple', marker='o')
plt.title('Median Values of Particulate Matter')
plt.xticks(x_ticks, timestamps, rotation=45)
plt.ylabel('Median PM2.5')
plt.legend()

# Standard deviation values
plt.subplot(3, 1, 3)
plt.plot(x_ticks, obs_stds, label='ObsValue Std Dev', color='blue', marker='o')
plt.plot(x_ticks, hofx0_stds, label='hofx0 Std Dev', color='green', marker='o')
plt.plot(x_ticks, hofx1_stds, label='hofx1 Std Dev', color='purple', marker='o')
plt.title('Standard Deviation of Particulate Matter')
plt.xticks(x_ticks, timestamps, rotation=45)
plt.ylabel('Std Dev PM2.5')
plt.legend()

plt.tight_layout()
plt.show()


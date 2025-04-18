import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import glob

# List of NetCDF files
file_pattern = './3dvar_diffusion_pm25_*.nc'
files = sorted(glob.glob(file_pattern))

# Variables to store the data
obs_values = []
hofx0_values = []
hofx1_values = []
timestamps = []

# Process each file
for file in files:
    with Dataset(file, 'r') as nc:
        # Extract time from the filename (assuming it's in the filename)
        timestamp = file.split('_')[-1].split('.')[0]  # e.g., 2020090118
        timestamps.append(timestamp)

        # Replace these with the correct variable names
        obs_var_name = 'ObsValue/particulatematter2p5Surface'  # Update this
        hofx0_var_name = 'hofx0/particulatematter2p5Surface'  # Update this
        hofx1_var_name = 'hofx1/particulatematter2p5Surface'  # Update this

        # Check for variable existence
        if obs_var_name in nc.variables:
            obs_values.append(nc.variables[obs_var_name][:])
        if hofx0_var_name in nc.variables:
            hofx0_values.append(nc.variables[hofx0_var_name][:])
        if hofx1_var_name in nc.variables:
            hofx1_values.append(nc.variables[hofx1_var_name][:])

# Convert to numpy arrays for easier manipulation
obs_values = np.concatenate(obs_values)
hofx0_values = np.concatenate(hofx0_values)
hofx1_values = np.concatenate(hofx1_values)

# Calculate statistics
obs_mean = np.mean(obs_values, axis=0)
obs_median = np.median(obs_values, axis=0)
obs_std = np.std(obs_values, axis=0)

hofx0_mean = np.mean(hofx0_values, axis=0)
hofx0_median = np.median(hofx0_values, axis=0)
hofx0_std = np.std(hofx0_values, axis=0)

hofx1_mean = np.mean(hofx1_values, axis=0)
hofx1_median = np.median(hofx1_values, axis=0)
hofx1_std = np.std(hofx1_values, axis=0)

# Create a figure for means
plt.figure(figsize=(15, 10))
plt.subplot(3, 1, 1)
plt.plot(timestamps, obs_mean, label='ObsValue Mean', color='blue')
plt.plot(timestamps, hofx0_mean, label='hofx0 Mean', color='green')
plt.plot(timestamps, hofx1_mean, label='hofx1 Mean', color='purple')
plt.title('Mean Values of Particulate Matter')
plt.xlabel('Time')
plt.ylabel('Mean PM2.5')
plt.xticks(rotation=45)
plt.legend()

# Create a figure for medians
plt.subplot(3, 1, 2)
plt.plot(timestamps, obs_median, label='ObsValue Median', color='blue')
plt.plot(timestamps, hofx0_median, label='hofx0 Median', color='green')
plt.plot(timestamps, hofx1_median, label='hofx1 Median', color='purple')
plt.title('Median Values of Particulate Matter')
plt.xlabel('Time')
plt.ylabel('Median PM2.5')
plt.xticks(rotation=45)
plt.legend()

# Create a figure for standard deviations
plt.subplot(3, 1, 3)
plt.plot(timestamps, obs_std, label='ObsValue Std Dev', color='blue')
plt.plot(timestamps, hofx0_std, label='hofx0 Std Dev', color='green')
plt.plot(timestamps, hofx1_std, label='hofx1 Std Dev', color='purple')
plt.title('Standard Deviation of Particulate Matter')
plt.xlabel('Time')
plt.ylabel('Std Dev PM2.5')
plt.xticks(rotation=45)
plt.legend()

plt.tight_layout()
plt.show()


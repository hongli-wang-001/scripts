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

        # Extract the variables
        obs_values.append(nc.variables['ObsValue/particulatematter2p5Surface'][:])
        hofx0_values.append(nc.variables['hofx0/particulatematter2p5Surface'][:])
        hofx1_values.append(nc.variables['hofx1/particulatematter2p5Surface'][:])

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

# Create time series plots
plt.figure(figsize=(15, 10))

# ObsValue
plt.subplot(3, 1, 1)
plt.plot(timestamps, obs_mean, label='Mean', color='blue')
plt.plot(timestamps, obs_median, label='Median', color='orange')
plt.fill_between(timestamps, obs_mean - obs_std, obs_mean + obs_std, color='blue', alpha=0.2)
plt.title('ObsValue/particulatematter2p5Surface')
plt.legend()

# Hofx0
plt.subplot(3, 1, 2)
plt.plot(timestamps, hofx0_mean, label='Mean', color='green')
plt.plot(timestamps, hofx0_median, label='Median', color='red')
plt.fill_between(timestamps, hofx0_mean - hofx0_std, hofx0_mean + hofx0_std, color='green', alpha=0.2)
plt.title('hofx0/particulatematter2p5Surface')
plt.legend()

# Hofx1
plt.subplot(3, 1, 3)
plt.plot(timestamps, hofx1_mean, label='Mean', color='purple')
plt.plot(timestamps, hofx1_median, label='Median', color='pink')
plt.fill_between(timestamps, hofx1_mean - hofx1_std, hofx1_mean + hofx1_std, color='purple', alpha=0.2)
plt.title('hofx1/particulatematter2p5Surface')
plt.legend()

plt.tight_layout()
plt.show()


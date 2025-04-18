import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Define file and variable details
netcdf_file = 'file1.nc'
groups = ['ObsValue', 'hofx0', 'hofx1']  # Replace with actual group names
variable_name = 'particulatematter2p5Surface'  # Actual variable name

# Function to extract data
def extract_data(nc_file, group_name, var_name):
    with Dataset(nc_file, 'r') as nc:
        data = nc.groups[group_name][var_name][:]
    return data

# Extract data for each group
data_dict = {}
for group in groups:
    data = extract_data(netcdf_file, group, variable_name)
    data_dict[group] = data

# Define bin range and number of bins
bin_range = (0, 120)
num_bins = 40
bins = np.linspace(bin_range[0], bin_range[1], num_bins + 1)

# Plot frequency distribution for original data
fig, axes = plt.subplots(nrows=1, ncols=len(groups), figsize=(15, 5), sharex='col', sharey='col')

for i, group in enumerate(groups):
    axes[i].hist(data_dict[group].flatten(), bins=bins, alpha=0.7, edgecolor='black')
    axes[i].set_title(f'{group}')
    axes[i].set_xlabel('Value')
    axes[i].set_ylabel('Frequency')
    axes[i].set_xlim(bin_range)

plt.tight_layout()
plt.savefig('frequency_distribution_original.png')
plt.show()


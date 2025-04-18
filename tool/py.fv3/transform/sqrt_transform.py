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

# Define bins for histogram
bin_size = 50
max_value = max(np.max(data) for data in data_dict.values())
bins = np.linspace(0, max_value, bin_size)

# Apply square root transform
transformed_data_dict = {group: np.sqrt(data) for group, data in data_dict.items()}

# Plot histograms before and after transformation
fig, axes = plt.subplots(nrows=2, ncols=len(groups), figsize=(15, 10), sharex='col', sharey='col')

for i, group in enumerate(groups):
    # Original Data
    axes[0, i].hist(data_dict[group].flatten(), bins=bins, alpha=0.7, edgecolor='black')
    axes[0, i].set_title(f'Original - {group}')
    axes[0, i].set_xlabel('Value')
    axes[0, i].set_ylabel('Frequency')

    # Transformed Data
    axes[1, i].hist(transformed_data_dict[group].flatten(), bins=bins, alpha=0.7, edgecolor='black')
    axes[1, i].set_title(f'Sqrt Transformed - {group}')
    axes[1, i].set_xlabel('Value')
    axes[1, i].set_ylabel('Frequency')

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig('frequency_distribution_before_after_sqrt_transform.png')

# Optional: Show plots
# plt.show()


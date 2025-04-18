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

# Apply log transform safely
transformed_data_dict = {}
epsilon = 1e-10  # Small constant to avoid log(0)
for group, data in data_dict.items():
    # Ensure no negative values
    data_non_negative = np.clip(data, epsilon, None)
    # Apply log transformation
    transformed_data = np.log(data_non_negative)
    transformed_data_dict[group] = transformed_data

# Define bins for original data
original_bin_range = (0, 120)
num_bins = 40
bins_original = np.linspace(original_bin_range[0], original_bin_range[1], num_bins + 1)  # +1 because np.linspace includes the endpoint

# Determine bins for transformed data
transformed_bin_ranges = {}
for group, transformed_data in transformed_data_dict.items():
    min_val = np.min(transformed_data)
    max_val = np.max(transformed_data)
    bins_transformed = np.linspace(min_val, max_val, num_bins + 1)
    transformed_bin_ranges[group] = bins_transformed

# Plot histograms before and after transformation
fig, axes = plt.subplots(nrows=2, ncols=len(groups), figsize=(15, 10), sharex='col', sharey='col')

for i, group in enumerate(groups):
    # Original Data
    axes[0, i].hist(data_dict[group].flatten(), bins=bins_original, alpha=0.7, edgecolor='black')
    axes[0, i].set_title(f'Original - {group}')
    axes[0, i].set_xlabel('Value')
    axes[0, i].set_ylabel('Frequency')
    axes[0, i].set_xlim(original_bin_range)

    # Transformed Data
    bins_transformed = transformed_bin_ranges[group]
    axes[1, i].hist(transformed_data_dict[group].flatten(), bins=bins_transformed, alpha=0.7, edgecolor='black')
    axes[1, i].set_title(f'Log Transformed - {group}')
    axes[1, i].set_xlabel('Log(Value)')
    axes[1, i].set_ylabel('Frequency')

# Adjust layout and save the figure
plt.tight_layout()
plt.savefig('frequency_distribution_before_after_log_transform_adjusted.png')

# Optional: Show plots
# plt.show()


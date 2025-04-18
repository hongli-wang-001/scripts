import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Define file and variable details
netcdf_file = 'file1.nc'
groups = ['ObsValue']  # Replace with actual group names
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

# Apply log transformation (handling zeros by adding a small constant)
transformed_data_dict = {}
for group, data in data_dict.items():
    # Replace zeros with a small positive number before log transformation
    data_cleaned = np.where(data > 0, data, np.nan)
    transformed_data_dict[group] = np.log1p(data_cleaned)

# Determine the bin range for log-transformed data
all_log_values = np.concatenate([transformed_data_dict[group].flatten() for group in groups])
all_log_values = all_log_values[~np.isnan(all_log_values)]  # Remove NaNs for bin calculation
bin_range_log = (0, np.max(all_log_values))
num_bins_log = 40
bins_log = np.linspace(bin_range_log[0], bin_range_log[1], num_bins_log + 1)

# Plot frequency distribution for log-transformed data
fig, axes = plt.subplots(nrows=1, ncols=len(groups), figsize=(15, 5), sharex='col', sharey='col')

for i, group in enumerate(groups):
    axes[i].hist(transformed_data_dict[group].flatten(), bins=bins_log, alpha=0.7, edgecolor='black')
    axes[i].set_title(f'{group} (Log Transformed)')
    axes[i].set_xlabel('Log(Value)')
    axes[i].set_ylabel('Frequency')
    axes[i].set_xlim(bin_range_log)

plt.tight_layout()
plt.savefig('frequency_distribution_log_transformed.png')
plt.show()


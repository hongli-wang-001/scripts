import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Define file and variable details
netcdf_file = 'file1.nc'  # Path to your NetCDF file
groups = ['ObsValue']  # List of groups in the NetCDF file
variable_name = 'particulatematter2p5Surface'  # Variable name to extract

# Function to extract data from NetCDF
def extract_data(nc_file, group_name, var_name):
    """Extract data from a specified NetCDF file, group, and variable."""
    with Dataset(nc_file, 'r') as nc:
        data = nc.groups[group_name][var_name][:]
    return data

# Extract data for each group
data_dict = {}
for group in groups:
    data = extract_data(netcdf_file, group, variable_name)
    data_dict[group] = data

# Apply log transformation, handling zeros by adding a small constant
transformed_data_dict = {}
for group, data in data_dict.items():
    # Replace zeros with NaN to avoid log(0) issues
    data_cleaned = np.where(data > 0, data, np.nan)
    # Apply log1p for a more numerically stable transformation (log(1 + x))
    transformed_data_dict[group] = np.log1p(data_cleaned)

# Determine the bin range for log-transformed data
all_log_values = np.concatenate([transformed_data_dict[group].flatten() for group in groups])
all_log_values = all_log_values[~np.isnan(all_log_values)]  # Remove NaNs for bin calculation
bin_range_log = (0, np.max(all_log_values))
num_bins_log = 40
bins_log = np.linspace(bin_range_log[0], bin_range_log[1], num_bins_log + 1)

# Compute original bin edges
bins_original = np.expm1(bins_log)  # Inverse of log1p

# Plot frequency distribution for log-transformed data
fig, axes = plt.subplots(nrows=1, ncols=len(groups), figsize=(15, 5), sharex='col', sharey='col')

# If only one subplot, axes will be a single Axes object
if len(groups) == 1:
    axes = [axes]

for i, group in enumerate(groups):
    ax = axes[i]
    # Plot log-transformed data
    counts, _ = np.histogram(transformed_data_dict[group].flatten(), bins=bins_log)
    ax.bar(bins_log[:-1], counts, width=np.diff(bins_log), alpha=0.7, edgecolor='black')
    
    # Setting primary x-axis (log-transformed values)
    ax.set_xlabel('Log(Value)')
    ax.set_ylabel('Frequency')
    ax.set_xlim(bin_range_log)
    
    # Adding secondary x-axis for original values
    ax2 = ax.twiny()
    ax2.set_xlim(np.log1p(bin_range_log[0]), np.log1p(bin_range_log[1]))
    ax2.set_xticks(np.log1p(bins_original))
    ax2.set_xticklabels([f'{int(val):,}' for val in bins_original])
    ax2.set_xlabel('Original Value')
    
    ax.set_title(f'{group} (Log Transformed)')
    
plt.tight_layout()
plt.savefig('frequency_distribution_log_transformed.png')
plt.show()


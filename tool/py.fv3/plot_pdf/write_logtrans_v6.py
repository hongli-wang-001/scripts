import numpy as np
from netCDF4 import Dataset

# Define file and variable details
netcdf_file = 'file1.nc'
new_netcdf_file = 'file_transformed_with_locations.nc'
groups = ['ObsValue', 'hofx0', 'hofx1']  # Replace with actual group names
variable_name = 'particulatematter2p5Surface'  # Actual variable name

# Define bin settings for log-transformed data
log_bin_range = (-2, 5)  # Bin range for log-transformed data
log_bin_size = 40  # Number of bins

# Function to extract data
def extract_data(nc_file, group_name, var_name):
    with Dataset(nc_file, 'r') as nc:
        data = nc.groups[group_name][var_name][:]
    return data

# Function to bin data and create bin indices
def create_bin_indices(data, bins):
    bin_indices = np.digitize(data.flatten(), bins) - 1
    bin_indices = np.clip(bin_indices, 0, len(bins) - 2)  # Ensure indices are within valid range
    return bin_indices.reshape(data.shape)

# Extract data for each group
data_dict = {}
for group in groups:
    data = extract_data(netcdf_file, group, variable_name)
    data_dict[group] = data

# Apply log transformation (handling zeros by adding a small constant)
transformed_data_dict = {}
for group, data in data_dict.items():
    data_cleaned = np.where(data > 0, data, np.nan)
    transformed_data_dict[group] = np.log1p(data_cleaned)

# Create bin edges for log-transformed data
log_min, log_max = log_bin_range
log_bins = np.linspace(log_min, log_max, log_bin_size + 1)
print(f"Log-transformed data bins: {log_bins}")

# Create a new NetCDF file to save the original and transformed data, and locations
with Dataset(new_netcdf_file, 'w', format='NETCDF4') as nc:
    for group in groups:
        data = data_dict[group]
        transformed_data = transformed_data_dict[group]

        data_shape = data.shape
        print(f"Creating dimensions for group {group}: {data_shape}")

        # Create dimensions based on the data shape
        for i, dim_size in enumerate(data_shape):
            nc.createDimension(f'dim_{group}_{i}', dim_size)

        # Create variable for original data
        var_name = f'{group}_{variable_name}'
        nc.createVariable(var_name, 'f4', [f'dim_{group}_{i}' for i in range(len(data_shape))])
        nc.variables[var_name][:] = data

        # Create variable for log-transformed data
        log_var_name = f'{group}_{variable_name}_log'
        nc.createVariable(log_var_name, 'f4', [f'dim_{group}_{i}' for i in range(len(data_shape))])
        nc.variables[log_var_name][:] = transformed_data

        # Create arrays to store bin indices for original and log-transformed data
        original_bin_indices = create_bin_indices(data, log_bins)
        log_bin_indices = create_bin_indices(transformed_data, log_bins)

        # Save bin indices to the NetCDF file
        for i in range(log_bin_size):
            # Original data bin indices
            original_bin_indices_var_name = f'{group}_{variable_name}_original_bin_{i + 1}'
            original_bin_var = nc.createVariable(original_bin_indices_var_name, 'i4', [f'dim_{group}_{d}' for d in range(len(data_shape))])
            original_bin_var[:] = (original_bin_indices == i).astype(int)

            # Log-transformed data bin indices
            log_bin_indices_var_name = f'{group}_{variable_name}_log_bin_{i + 1}'
            log_bin_var = nc.createVariable(log_bin_indices_var_name, 'i4', [f'dim_{group}_{d}' for d in range(len(data_shape))])
            log_bin_var[:] = (log_bin_indices == i).astype(int)

        print(f"Original data, log-transformed data, and bin indices have been saved to {new_netcdf_file}.")


import numpy as np
from netCDF4 import Dataset

# Define file and variable details
netcdf_file = 'file1.nc'
new_netcdf_file = 'file_transformed.nc'
groups = ['ObsValue', 'hofx0', 'hofx1']  # Replace with actual group names
variable_name = 'particulatematter2p5Surface'  # Actual variable name

# Define bins for log-transformed data
bin_range_log = (-2, 5)  # Adjust as needed
bin_size_log = 40
bins_log = np.linspace(bin_range_log[0], bin_range_log[1], bin_size_log)

# Function to extract data
def extract_data(nc_file, group_name, var_name):
    with Dataset(nc_file, 'r') as nc:
        data = nc.groups[group_name][var_name][:]
    return data

# Function to calculate bin locations
def calculate_bin_locations(data, bins):
    bin_indices = np.digitize(data.ravel(), bins) - 1
    bin_locations = np.full(data.shape, np.nan)
    for i in range(len(bins) - 1):
        bin_locations.ravel()[bin_indices == i] = i
    return bin_locations

# Extract data for each group
data_dict = {}
for group in groups:
    data = extract_data(netcdf_file, group, variable_name)
    data_dict[group] = data

# Apply log transformation (handling zeros by adding a small constant)
transformed_data_dict = {}
bin_locations_original_dict = {}
bin_locations_transformed_dict = {}
for group, data in data_dict.items():
    # Replace zeros with a small positive number before log transformation
    data_cleaned = np.where(data > 0, data, np.nan)
    transformed_data = np.log1p(data_cleaned)
    
    # Calculate bin locations for original and transformed data
    bin_locations_original = calculate_bin_locations(data, bins_log)
    bin_locations_transformed = calculate_bin_locations(transformed_data, bins_log)
    
    transformed_data_dict[group] = transformed_data
    bin_locations_original_dict[group] = bin_locations_original
    bin_locations_transformed_dict[group] = bin_locations_transformed

# Create a new NetCDF file to save the original, transformed data, and bin locations
with Dataset(new_netcdf_file, 'w', format='NETCDF4') as nc:
    for group in groups:
        data = data_dict[group]
        data_shape = data.shape
        
        # Create dimensions dynamically
        for dim, size in enumerate(data_shape):
            nc.createDimension(f'dim_{group}_{dim}', size)

        # Create variable for original data
        var_name = f'{group}_{variable_name}'
        nc.createVariable(var_name, 'f4', [f'dim_{group}_{dim}' for dim in range(len(data_shape))])
        nc.variables[var_name][:] = data

        # Create variable for log-transformed data
        log_var_name = f'{group}_{variable_name}_log'
        nc.createVariable(log_var_name, 'f4', [f'dim_{group}_{dim}' for dim in range(len(data_shape))])
        nc.variables[log_var_name][:] = transformed_data_dict[group]
        
        # Create variable for bin locations in original data
        bin_loc_orig_name = f'{group}_{variable_name}_bin_loc'
        nc.createVariable(bin_loc_orig_name, 'i4', [f'dim_{group}_{dim}' for dim in range(len(data_shape))])
        nc.variables[bin_loc_orig_name][:] = bin_locations_original_dict[group]

        # Create variable for bin locations in log-transformed data
        bin_loc_log_name = f'{group}_{variable_name}_log_bin_loc'
        nc.createVariable(bin_loc_log_name, 'i4', [f'dim_{group}_{dim}' for dim in range(len(data_shape))])
        nc.variables[bin_loc_log_name][:] = bin_locations_transformed_dict[group]

    print(f"Original, log-transformed data, and bin location information have been saved to {new_netcdf_file}.")


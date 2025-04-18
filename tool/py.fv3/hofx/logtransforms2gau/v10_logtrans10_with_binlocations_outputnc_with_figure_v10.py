import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt

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

# Function to plot histogram
def plot_histogram(data, bins, title, filename):
    plt.figure(figsize=(10, 6))
    plt.hist(data.ravel(), bins=bins, edgecolor='k', alpha=0.7)
    plt.title(title)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.xlim(min(bins), max(bins))
    plt.grid(True)
    plt.savefig(filename)
    plt.close()

# Extract data for each group
data_dict = {}
for group in groups:
    data = extract_data(netcdf_file, group, variable_name)
    data_dict[group] = data

# Apply log transformation (adding 1 to handle zeros)
transformed_data_dict = {}
bin_locations_original_dict = {}
bin_locations_transformed_dict = {}
for group, data in data_dict.items():
    # Add 1 before applying log transformation
    data_cleaned = np.where(data >= 0, data + 1, 1)  # Add 1 to handle zeros
    transformed_data = np.log10(data_cleaned)  # Base-10 log transformation
    
    # Calculate bin locations for original and transformed data
    bin_locations_original = calculate_bin_locations(data, bins_log)
    bin_locations_transformed = calculate_bin_locations(transformed_data, bins_log)
    
    transformed_data_dict[group] = transformed_data
    bin_locations_original_dict[group] = bin_locations_original
    bin_locations_transformed_dict[group] = bin_locations_transformed

    # Plot the frequency of the log-transformed data
    plot_histogram(transformed_data, bins_log, f'Log-Transformed Data - {group}', f'{group}_{variable_name}_log_histogram.png')

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


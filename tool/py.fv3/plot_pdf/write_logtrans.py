import numpy as np
from netCDF4 import Dataset

# Define file and variable details
netcdf_file = 'file1.nc'
new_netcdf_file = 'file_transformed.nc'
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

# Apply log transformation (handling zeros by adding a small constant)
transformed_data_dict = {}
for group, data in data_dict.items():
    # Replace zeros with a small positive number before log transformation
    data_cleaned = np.where(data > 0, data, np.nan)
    transformed_data_dict[group] = np.log1p(data_cleaned)

# Create a new NetCDF file to save the original and transformed data
with Dataset(new_netcdf_file, 'w', format='NETCDF4') as nc:
    # Create dimensions
    for group in groups:
        data_shape = data_dict[group].shape
        nc.createDimension('dim_' + group, data_shape[0])
        nc.createDimension('dim_y', data_shape[1])
        nc.createDimension('dim_x', data_shape[2])

    # Create variables for original and transformed data
    for group in groups:
        # Create variable for original data
        nc.createVariable(f'{group}_{variable_name}', 'f4', ('dim_' + group, 'dim_y', 'dim_x'))
        nc.variables[f'{group}_{variable_name}'][:] = data_dict[group]

        # Create variable for log-transformed data
        nc.createVariable(f'{group}_{variable_name}_log', 'f4', ('dim_' + group, 'dim_y', 'dim_x'))
        nc.variables[f'{group}_{variable_name}_log'][:] = transformed_data_dict[group]

    print(f"Original and log-transformed data have been saved to {new_netcdf_file}.")


import netCDF4 as nc
import numpy as np
from scipy.ndimage import gaussian_filter
import os

# Define the input and output NetCDF file names and the variable name
input_netcdf_file = 'file1_fcst.nc'
output_netcdf_file = 'filtered_output.nc'
variable_name = 'asvoo1j'

# Define the Gaussian filter standard deviation
sigma = 5.0  # Adjust this value as needed for your filtering

# Function to check file existence and print path
def check_file_existence(file_path):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    print("File path:", file_path)

# Check if the input file exists
check_file_existence(input_netcdf_file)

# Open the input NetCDF file and read the variable
with nc.Dataset(input_netcdf_file, 'r') as ds:
    # Check if the variable exists in the file
    if variable_name not in ds.variables:
        raise KeyError(f"The variable '{variable_name}' is not present in the file.")
    
    # Read the 4D array from the specified variable
    data = ds.variables[variable_name][:]
    
    # Ensure the data is 4D
    if data.ndim != 4:
        raise ValueError("The variable does not contain a 4D array.")
    
    # Get dimensions
    time_dim, z_dim, y_dim, x_dim = data.shape

    # Initialize an array for the filtered data
    filtered_data = np.empty_like(data)

    # Apply the Gaussian filter to each 2D slice in the time dimension
    for t in range(time_dim):
        for z in range(z_dim):
            filtered_data[t, z, :, :] = gaussian_filter(data[t, z, :, :], sigma=sigma)

    # Save the filtered data to a new NetCDF file
    with nc.Dataset(output_netcdf_file, 'w', format='NETCDF4') as ds_out:
        # Create dimensions
        ds_out.createDimension('time', time_dim)
        ds_out.createDimension('zaxis_1', z_dim)
        ds_out.createDimension('yaxis_1', y_dim)
        ds_out.createDimension('xaxis_1', x_dim)

        # Create variable
        filtered_var = ds_out.createVariable(variable_name, 'f4', ('time', 'zaxis_1', 'yaxis_1', 'xaxis_1'))

        # Write the filtered data
        filtered_var[:] = filtered_data

        # Copy attributes if needed
        try:
            for attr in ds.variables[variable_name].ncattrs():
                filtered_var.setncattr(attr, ds.variables[variable_name].getncattr(attr))
        except Exception as e:
            print(f"Error copying attributes: {e}")

print(f"Filtered 4D data saved to {output_netcdf_file}")


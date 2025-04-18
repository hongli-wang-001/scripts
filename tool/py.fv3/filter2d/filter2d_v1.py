import netCDF4 as nc
import numpy as np
import scipy.ndimage as ndimage
from scipy.ndimage import gaussian_filter
from netCDF4 import Dataset

# Define the input and output NetCDF file names and the variable name
input_netcdf_file = 'input_file.nc'
output_netcdf_file = 'filtered_output.nc'
variable_name = 'variable'

# Define the Gaussian filter standard deviation
sigma = 2.0  # Adjust this value as needed for your filtering

# Open the input NetCDF file and read the variable
with nc.Dataset(input_netcdf_file, 'r') as ds:
    # Read the 2D array from the specified variable
    data = ds.variables[variable_name][:]
    
    # Ensure the data is 2D
    if data.ndim != 2:
        raise ValueError("The variable does not contain a 2D array.")

# Apply a Gaussian filter to the 2D array
filtered_data = gaussian_filter(data, sigma=sigma)

# Save the filtered data to a new NetCDF file
with nc.Dataset(output_netcdf_file, 'w', format='NETCDF4') as ds_out:
    # Create dimensions
    ds_out.createDimension('x', data.shape[0])
    ds_out.createDimension('y', data.shape[1])
    
    # Create variable
    filtered_var = ds_out.createVariable(variable_name, 'f4', ('x', 'y'))
    
    # Write the filtered data
    filtered_var[:] = filtered_data
    
    # Copy attributes if needed
    for attr in ds.variables[variable_name].ncattrs():
        filtered_var.setncattr(attr, ds.variables[variable_name].getncattr(attr))

print(f"Filtered data saved to {output_netcdf_file}")


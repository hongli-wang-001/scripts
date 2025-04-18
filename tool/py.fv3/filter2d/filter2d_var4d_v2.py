import netCDF4 as nc
import numpy as np
import scipy.ndimage as ndimage
from scipy.ndimage import gaussian_filter

# Define the input and output NetCDF file names and the variable name
input_netcdf_file = 'file1.nc'
output_netcdf_file = 'filtered_output.nc'
variable_name = 'asoil'

# Define the Gaussian filter standard deviation
sigma = 2.0  # Adjust this value as needed for your filtering

# Open the input NetCDF file and read the variable
with nc.Dataset(input_netcdf_file, 'r') as ds:
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
    for attr in ds.variables[variable_name].ncattrs():
        filtered_var.setncattr(attr, ds.variables[variable_name].getncattr(attr))

print(f"Filtered 4D data saved to {output_netcdf_file}")


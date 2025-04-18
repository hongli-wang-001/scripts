import os
import netCDF4 as nc
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.ndimage import median_filter
from scipy.ndimage import uniform_filter

# Define the input and output NetCDF file names and the list of variables
input_netcdf_file = 'file1_fcst.nc'
output_netcdf_file = 'filtered_mean_6by6_output.nc'
variables = [
    'aso4i', 'ano3i', 'anh4i', 'anai', 'acli', 'aeci', 'aothri', 'alvpo1i',
    'asvpo1i', 'asvpo2i', 'alvoo1i', 'alvoo2i', 'asvoo1i', 'asvoo2i', 'aso4j',
    'ano3j', 'anh4j', 'anaj', 'aclj', 'aecj', 'aothrj', 'afej', 'asij', 'atij',
    'acaj', 'amgj', 'amnj', 'aalj', 'akj', 'alvpo1j', 'asvpo1j', 'asvpo2j',
    'asvpo3j', 'aivpo1j', 'axyl1j', 'axyl2j', 'axyl3j', 'atol1j', 'atol2j',
    'atol3j', 'abnz1j', 'abnz2j', 'abnz3j', 'aiso1j', 'aiso2j', 'aiso3j',
    'atrp1j', 'atrp2j', 'asqtj', 'aalk1j', 'aalk2j', 'apah1j', 'apah2j',
    'apah3j', 'aorgcj', 'aolgbj', 'aolgaj', 'alvoo1j', 'alvoo2j', 'asvoo1j',
    'asvoo2j', 'asvoo3j', 'apcsoj', 'aso4k', 'asoil', 'acors', 'aseacat',
    'aclk', 'ano3k', 'anh4k'
]

# Define the Gaussian filter standard deviation
sigma = 3.0  # Adjust this value as needed for your filtering

# Function to apply log transformation
def log_transform(data):
    """Apply log transformation to data."""
    return np.log1p(data)  # Using np.log1p for log(1 + data) to handle zero values

# Function to apply inverse log transformation
def inverse_log_transform(data):
    """Apply inverse log transformation to data."""
    return np.expm1(data)  # Using np.expm1 for exp(data) - 1

# Function to check file existence and print path
def check_file_existence(file_path):
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file {file_path} does not exist.")
    print("File path:", file_path)

def process_variable(input_netcdf_file, variable_name, sigma, output_file):
    try:
        # Open the input NetCDF file and read the variable
        with nc.Dataset(input_netcdf_file, 'r') as ds:
            # Check if the variable exists
            if variable_name not in ds.variables:
                raise ValueError(f"Variable '{variable_name}' not found in the input file.")

            # Read the 4D array from the specified variable
            data = ds.variables[variable_name][:]

            # Ensure the data is 4D
            if data.ndim != 4:
                raise ValueError(f"The variable '{variable_name}' does not contain a 4D array.")

            # Get dimensions
            time_dim, z_dim, y_dim, x_dim = data.shape

        # Initialize an array for the filtered data
        filtered_data = np.empty_like(data)

        # Apply the log transformation, Gaussian filter, and inverse log transformation
        #for t in range(time_dim):
        #    for z in range(z_dim):
        #        # Apply log transformation
        #        log_data = log_transform(data[t, z, :, :])

                # Apply Gaussian filter
        #        filtered_log_data = gaussian_filter(log_data, sigma=sigma)

                # Apply inverse log transformation
        #        filtered_data[t, z, :, :] = inverse_log_transform(filtered_log_data)

#filtered_data = median_filter(data, size=(3, 3))
        # Apply the median filter to each 2D slice in the time and vertical dimensions
        for t in range(time_dim):
            for z in range(z_dim):
                filtered_data[t, z, :, :] = uniform_filter(data[t, z, :, :], size=(6, 6))

        # Save the filtered data to a new NetCDF file
        mode = 'a' if os.path.isfile(output_file) else 'w'
        with nc.Dataset(output_file, mode, format='NETCDF4') as ds_out:
            # Create dimensions if they do not exist
            if 'time' not in ds_out.dimensions:
                ds_out.createDimension('time', time_dim)
            if 'zaxis_1' not in ds_out.dimensions:
                ds_out.createDimension('zaxis_1', z_dim)
            if 'yaxis_1' not in ds_out.dimensions:
                ds_out.createDimension('yaxis_1', y_dim)
            if 'xaxis_1' not in ds_out.dimensions:
                ds_out.createDimension('xaxis_1', x_dim)

            # Create or update variable
            if variable_name in ds_out.variables:
                filtered_var = ds_out.variables[variable_name]
            else:
                filtered_var = ds_out.createVariable(variable_name, 'f4', ('time', 'zaxis_1', 'yaxis_1', 'xaxis_1'))

            # Write the filtered data
            filtered_var[:] = filtered_data

            # Copy attributes from the original variable
            with nc.Dataset(input_netcdf_file, 'r') as ds:
                for attr in ds.variables[variable_name].ncattrs():
                    filtered_var.setncattr(attr, ds.variables[variable_name].getncattr(attr))

        print(f"Filtered data for variable '{variable_name}' saved to {output_file}")

    except FileNotFoundError:
        print(f"Error: The file {input_netcdf_file} does not exist.")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Check if the input file exists
check_file_existence(input_netcdf_file)

# Process each variable
for variable in variables:
    process_variable(input_netcdf_file, variable, sigma, output_netcdf_file)


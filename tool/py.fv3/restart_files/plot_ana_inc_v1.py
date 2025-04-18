import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Define the input NetCDF file name and the list of variables
input_netcdf_file = 'inc_20200902.00.nc'
output_directory = 'figures'
print("Current working directory:", os.getcwd())
print("File path:", input_netcdf_file)
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

# Ensure output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

def plot_variable_at_levels(data, variable_name, z_dim, output_directory):
    """
    Plots the 2D slices of a 4D variable at different zaxis_1 levels.
    """
    time_dim, _, y_dim, x_dim = data.shape

    # Iterate over each zaxis_1 level
    for z in range(z_dim):
        plt.figure(figsize=(12, 6))
        plt.title(f'{variable_name} - zaxis_1 level {z}')
        plt.imshow(data[0, z, :, :], cmap='viridis', origin='lower')  # Assuming time_dim = 1
        plt.colorbar(label='Value')
        plt.xlabel('xaxis_1')
        plt.ylabel('yaxis_1')
        plt.savefig(os.path.join(output_directory, f'{variable_name}_z{z}.png'))
        plt.close()

def main():
    try:
        # Open the input NetCDF file and read the variables
        with nc.Dataset(input_netcdf_file, 'r') as ds:
            for variable_name in variables:
                if variable_name in ds.variables:
                    data = ds.variables[variable_name][:]
                    if data.ndim == 4:
                        time_dim, z_dim, y_dim, x_dim = data.shape
                        # Plot data for each zaxis_1 level
                        plot_variable_at_levels(data, variable_name, z_dim, output_directory)
                    else:
                        print(f"Variable '{variable_name}' does not have 4 dimensions.")
                else:
                    print(f"Variable '{variable_name}' not found in the input file.")
    except FileNotFoundError:
        print(f"Error: The file {input_netcdf_file} does not exist.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    main()


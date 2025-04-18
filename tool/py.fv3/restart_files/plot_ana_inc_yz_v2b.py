import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Define the input NetCDF file name and the list of variables
input_netcdf_file = 'inc_20200902.00.nc'
output_directory = 'figures'
xaxis_indices = [42, 57]  # Indices of xaxis_1 to plot
print("Current working directory:", os.getcwd())
print("File path:", input_netcdf_file)

variables = ['aso4j','ano3j', 'anh4j','aecj','aorgcj','asoil', 'acors', 'aseacat']

# Ensure output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

def plot_slices_at_xaxis(data, variable_name, xaxis_indices, output_directory):
    """
    Plots the 2D slices of a 4D variable at specified xaxis_1 values.
    """
    time_dim, z_dim, y_dim, x_dim = data.shape

    # Iterate over each specified xaxis_1 index
    for x in xaxis_indices:
        if x < x_dim:
            # Extract 2D slices for all zaxis_1 levels at the specified xaxis_1 index
            slice_data = data[0, :, :, x]  # Assuming time_dim = 1

            # Check the shape of the slice
            if slice_data.ndim == 3:
                # Plot each zaxis_1 level
                for z in range(slice_data.shape[0]):
                    plt.figure(figsize=(12, 6))
                    plt.title(f'{variable_name} - xaxis_1 level {x} - zaxis_1 level {z}')
                    plt.imshow(slice_data[z, :, :], cmap='viridis', origin='lower')
                    plt.colorbar(label='Value')
                    plt.xlabel('yaxis_1')
                    plt.ylabel('zaxis_1')
                    plt.savefig(os.path.join(output_directory, f'{variable_name}_x{x}_z{z}.png'))
                    plt.close()
            else:
                print(f"Warning: Extracted slice for xaxis_1 = {x} has unexpected shape {slice_data.shape}.")
        else:
            print(f"Warning: xaxis_1 index {x} is out of bounds for variable '{variable_name}'.")

def main():
    try:
        # Open the input NetCDF file and read the variables
        with nc.Dataset(input_netcdf_file, 'r') as ds:
            for variable_name in variables:
                if variable_name in ds.variables:
                    data = ds.variables[variable_name][:]
                    if data.ndim == 4:
                        time_dim, z_dim, y_dim, x_dim = data.shape
                        # Plot data for each xaxis_1 index
                        plot_slices_at_xaxis(data, variable_name, xaxis_indices, output_directory)
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


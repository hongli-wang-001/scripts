import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Define the input NetCDF file name and the list of variables
input_netcdf_file = 'bkg.nc'
output_directory = 'figures_bkg_yz'
xaxis_indices = [42, 57]  # Indices of xaxis_1 to plot
print("Current working directory:", os.getcwd())
print("File path:", input_netcdf_file)
variables = ['aso4j','ano3j', 'anh4j','aecj','aorgcj','asoil', 'acors', 'aseacat']

# Ensure output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

def create_symmetric_cmap():
    """Create a symmetric colormap centered on white."""
    colors = ['white', 'blue', 'red']
    cmap = mcolors.LinearSegmentedColormap.from_list('symmetric_cmap', colors)
    return cmap

def plot_slices_at_xaxis(data, variable_name, xaxis_indices, output_directory):
    """
    Plots the 2D slices of a 4D variable at specified xaxis_1 values.
    Reverses the z-direction of the plot.
    """
    time_dim, z_dim, y_dim, x_dim = data.shape
    #cmap = create_symmetric_cmap()
    cmap='YlOrBr'
    # Iterate over each specified xaxis_1 index
    for x in xaxis_indices:
        if x < x_dim:
            # Extract the 2D slice for the specified xaxis_1 index
            slice_data = data[0, :, :, x]  # Assuming time_dim = 1
            
            # Reverse the z-direction (flip vertically)
            slice_data_reversed = np.flipud(slice_data)

            # Find the symmetric limits around zero
            min_val = np.min(slice_data_reversed)
            max_val = np.max(slice_data_reversed)
            #vmin = min(min_val, -max_val)
            #vmax = max(max_val, -min_val)
            vmin = 0.01
            vmax = max_val 
            # scale vmin and vmax
            #vmin = vmin/2.0
            vmax = vmax/2.0
 
            # Plot the 2D slice
            plt.figure(figsize=(12, 6))
            plt.title(f'{variable_name} - xaxis_1 level {x}')
            plt.imshow(slice_data_reversed, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)
            plt.colorbar(label='Value')
            plt.xlabel('yaxis_1')
            plt.ylabel('zaxis_1')
            plt.savefig(os.path.join(output_directory, f'{variable_name}_x{x}_reversed_symmetric.png'))
            plt.close()
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


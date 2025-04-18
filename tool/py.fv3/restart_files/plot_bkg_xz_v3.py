import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Define the input NetCDF file name and the list of variables
input_netcdf_file = 'bkg.nc'
output_directory = 'figures_bkg_xz'
yaxis_indices = [111, 147]  # Indices of yaxis_1 to plot
xaxis_indices = [42, 57]  # Indices of xaxis_1 to plot
print("Current working directory:", os.getcwd())
print("File path:", input_netcdf_file)
variables = ['aso4j', 'ano3j', 'anh4j', 'aecj', 'aorgcj', 'asoil', 'acors', 'aseacat']

# Ensure output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

def create_symmetric_cmap():
    """Create a symmetric colormap centered on white."""
    colors = ['white', 'blue', 'red']
    cmap = mcolors.LinearSegmentedColormap.from_list('symmetric_cmap', colors)
    return cmap

def plot_slices_at_yaxis(data, variable_name, yaxis_indices, output_directory):
    """
    Plots the 2D slices of a 4D variable at specified yaxis_1 values across all xaxis_1 values.
    Reverses the z-direction of the plot.
    """
    time_dim, z_dim, y_dim, x_dim = data.shape
    cmap = create_symmetric_cmap()  # Use the symmetric colormap

    # Iterate over each specified yaxis_1 index
    for y in yaxis_indices:
        if y < y_dim:
            # Extract the 2D slice for the specified yaxis_1 index
            slice_data = data[0, :, y, :]  # Assuming time_dim = 1

            # Reverse the z-direction (flip vertically)
            slice_data_reversed = np.flipud(slice_data)

            # Find the symmetric limits around zero
            min_val = np.min(slice_data_reversed)
            max_val = np.max(slice_data_reversed)
            vmin = min_val
            vmax = max_val

            # Plot the 2D slice
            plt.figure(figsize=(12, 6))
            plt.title(f'{variable_name} - yaxis_1 level {y}')
            img = plt.imshow(slice_data_reversed, cmap=cmap, origin='lower', vmin=vmin, vmax=vmax)

            # Add colorbar with adjusted size
            cbar = plt.colorbar(img, shrink=0.7)
            cbar.set_label('Value')

            # Print min, max, mean, median, and std on the figure
            stats_text = (f'Min: {min_val:.2f}\n'
                          f'Max: {max_val:.2f}\n'
                          f'Mean: {np.mean(slice_data_reversed):.2f}\n'
                          f'Median: {np.median(slice_data_reversed):.2f}\n'
                          f'Std: {np.std(slice_data_reversed):.2f}')
            plt.gca().text(0.05, 0.05, stats_text, transform=plt.gca().transAxes,
                           fontsize=12, verticalalignment='bottom', bbox=dict(boxstyle='round', alpha=0.5, facecolor='white'))

            plt.xlabel('xaxis_1')
            plt.ylabel('zaxis_1')
            plt.savefig(os.path.join(output_directory, f'{variable_name}_y{y}_reversed.png'))
            plt.close()
        else:
            print(f"Warning: yaxis_1 index {y} is out of bounds for variable '{variable_name}'.")

def main():
    try:
        # Open the input NetCDF file and read the variables
        with nc.Dataset(input_netcdf_file, 'r') as ds:
            for variable_name in variables:
                if variable_name in ds.variables:
                    data = ds.variables[variable_name][:]
                    if data.ndim == 4:
                        time_dim, z_dim, y_dim, x_dim = data.shape
                        # Plot data for each yaxis_1 index
                        plot_slices_at_yaxis(data, variable_name, yaxis_indices, output_directory)
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


import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define the input NetCDF file name and the list of variables
input_netcdf_file = 'inc_20200902.00.nc'
grid_file = 'fv3_grid_spec.nc'  # Assuming you have a separate file for grid lat/lon
output_directory = 'figures_xy'
print("Current working directory:", os.getcwd())
print("File path:", input_netcdf_file)
variables = ['aso4j', 'ano3j', 'anh4j', 'aecj', 'aorgcj', 'asoil', 'acors', 'aseacat']

# Ensure output directory exists
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

def create_symmetric_cmap():
    """Create a symmetric colormap centered on white."""
    import matplotlib.colors as mcolors
    colors = ['blue', 'white', 'red']
    cmap = mcolors.LinearSegmentedColormap.from_list('symmetric_cmap', colors)
    return cmap

def plot_variable_at_levels(data, variable_name, z_dim, lats, lons, output_directory):
    """
    Plots the 2D slices of a 4D variable at different zaxis_1 levels with latitude and longitude.
    """
    time_dim, _, y_dim, x_dim = data.shape
    cmap = create_symmetric_cmap()

    # Iterate over each zaxis_1 level
    for z in range(z_dim):
        plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        plt.title(f'{variable_name} - zaxis_1 level {z}')
        
        # Plot data
        img = ax.imshow(data[0, z, :, :], cmap=cmap, origin='lower',
                        vmin=-np.max(data[0, z, :, :]), vmax=np.max(data[0, z, :, :]),
                        extent=[np.min(lons), np.max(lons), np.min(lats), np.max(lats)],
                        transform=ccrs.PlateCarree())
        
        # Add map features
        ax.add_feature(cfeature.COASTLINE, linestyle='--')
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.STATES, linestyle='-')
        ax.gridlines(draw_labels=True)
        
        # Add colorbar
        plt.colorbar(img, ax=ax, label='Value')
        
        # Save the figure
        plt.savefig(os.path.join(output_directory, f'{variable_name}_z{z}_map.png'))
        plt.close()

def main():
    try:
        # Open the input NetCDF file and read the variables
        with nc.Dataset(input_netcdf_file, 'r') as ds:
            # Open the grid file to get lat/lon
            with nc.Dataset(grid_file, 'r') as grid_ds:
                lats = grid_ds.variables['grid_lat'][:]  # Replace with actual variable names
                lons = grid_ds.variables['grid_lon'][:]  # Replace with actual variable names
            
            for variable_name in variables:
                if variable_name in ds.variables:
                    data = ds.variables[variable_name][:]
                    if data.ndim == 4:
                        time_dim, z_dim, y_dim, x_dim = data.shape
                        # Plot data for each zaxis_1 level
                        plot_variable_at_levels(data, variable_name, z_dim, lats, lons, output_directory)
                    else:
                        print(f"Variable '{variable_name}' does not have 4 dimensions.")
                else:
                    print(f"Variable '{variable_name}' not found in the input file.")
    except FileNotFoundError:
        print(f"Error: The file {input_netcdf_file} or {grid_file} does not exist.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    main()


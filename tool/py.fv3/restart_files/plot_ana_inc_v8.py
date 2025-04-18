import os
import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Define the input NetCDF file name and the list of variables
input_netcdf_file = 'inc_20200902.00.nc'
grid_file = 'fv3_grid_spec.nc'  # File containing grid lat/lon
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

def wrap_longitudes(lons):
    """Wrap longitudes to the range [0, 360)."""
    return np.where(lons < 0, lons + 360, lons)

def print_statistics(data, z):
    """Print statistical information of the 2D data."""
    print(f'\nStatistics for zaxis_1 level {z}:')
    print(f'Min: {np.min(data)}')
    print(f'Max: {np.max(data)}')
    print(f'Mean: {np.mean(data)}')
    print(f'Median: {np.median(data)}')
    print(f'Standard Deviation: {np.std(data)}')

def plot_variable_at_levels(data, variable_name, z_dim, lats, lons, output_directory):
    """
    Plots the 2D slices of a 4D variable at every 5th zaxis_1 level with latitude and longitude.
    """
    time_dim, _, y_dim, x_dim = data.shape
    cmap = create_symmetric_cmap()

    # Define the extent for zooming into the United States
    extent = [-130, -60, 24, 50]  # [min_lon, max_lon, min_lat, max_lat]
    
    # Normalize longitudes if needed
    lons = wrap_longitudes(lons)
    
    # Iterate over every 5th zaxis_1 level
    for z in range(63, z_dim, 5):
        plt.figure(figsize=(12, 8))
        ax = plt.axes(projection=ccrs.PlateCarree())
        plt.title(f'{variable_name} - zaxis_1 level {z}')
        
        # Extract 2D data for plotting
        slice_data = data[0, z, :, :]
        
        # Print statistics
        print_statistics(slice_data, z)

        # Find the symmetric limits around zero
        min_val = np.min(slice_data)
        max_val = np.max(slice_data)
        vmin = min(min_val, -max_val)
        vmax = max(max_val, -min_val)

        # scale vmin and vmax
        vmin = vmin/5.0
        vmax = vmax/5.0
        vmin = max(vmin, -100)
        vmax = min(vmax, 100)

        # Compute standard deviation
        std_dev = np.std(slice_data)
        
        # Set vmin and vmax based on ±3 * std deviation
        vmin = -5 * std_dev
        vmax = 5 * std_dev
        vmin = max(vmin, -100)
        vmax = min(vmax, 100)
 
        # Plot data
        img = ax.imshow(slice_data, cmap=cmap, origin='lower',
#                        vmin=-np.max(slice_data), vmax=np.max(slice_data),
                        vmin=vmin, vmax=vmax,
                        extent=extent,
                        transform=ccrs.PlateCarree())
        
        # Add built-in map features without downloading new data
        ax.add_feature(cfeature.COASTLINE, linestyle='--')
        ax.add_feature(cfeature.BORDERS, linestyle=':')
        ax.add_feature(cfeature.STATES, linestyle='-')
        ax.gridlines(draw_labels=True)
        
        # Add colorbar
        # plt.colorbar(img, ax=ax, label='Value')
        cbar = plt.colorbar(img, ax=ax, orientation='vertical', fraction=0.018, pad=0.06)
        cbar.set_label('Value')

        # Add statistics text in the bottom-right corner
        stats_text = (f"Min: {np.min(slice_data):.2f}\n"
                      f"Max: {np.max(slice_data):.2f}\n"
                      f"Mean: {np.mean(slice_data):.2f}\n"
                      f"Median: {np.median(slice_data):.2f}\n"
                      f"Std: {np.std(slice_data):.2f}")
        
        plt.text(0.95, 0.05, stats_text, transform=ax.transAxes,
                 fontsize=10, verticalalignment='bottom', horizontalalignment='right',
                 bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
                
        # Save the figure
        plt.savefig(os.path.join(output_directory, f'{variable_name}_z{z}_map.png'))
        plt.close()

def main():
    try:
        # Open the input NetCDF file and read the variables
        with nc.Dataset(input_netcdf_file, 'r') as ds:
            # Open the grid file to get lat/lon
            with nc.Dataset(grid_file, 'r') as grid_ds:
                lats = grid_ds.variables['grid_latt'][:]  # 2D latitude array
                lons = grid_ds.variables['grid_lont'][:]  # 2D longitude array
            
            for variable_name in variables:
                if variable_name in ds.variables:
                    data = ds.variables[variable_name][:]
                    if data.ndim == 4:
                        time_dim, z_dim, y_dim, x_dim = data.shape
                        # Plot data for every 5th zaxis_1 level
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


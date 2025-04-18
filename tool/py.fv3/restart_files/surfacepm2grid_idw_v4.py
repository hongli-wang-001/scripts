import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors

def inverse_distance_weighting(xy_obs, values_obs, xy_grid, cutoff_radius=0.2, power=2):
    tree = cKDTree(xy_obs)
    distances, indices = tree.query(xy_grid, k=len(xy_obs), distance_upper_bound=cutoff_radius)
    
    interpolated_values = np.full(len(xy_grid), np.nan)  # Initialize with NaNs
    for i, (dist, idx) in enumerate(zip(distances, indices)):
        # Filter out infinite distances
        valid_idx = idx[dist < np.inf]
        valid_dist = dist[dist < np.inf]
        
        if valid_idx.size == 0:
            print(f"Warning: No valid observations within cutoff radius for grid point {i}.")
            continue
        
        # Calculate weights
        weights = 1 / (valid_dist ** power + np.finfo(float).eps)  # Avoid division by zero
        
        # Debugging prints
        print(f"Grid point {i}:")
        print(f"  Valid indices: {valid_idx}")
        print(f"  Valid distances: {valid_dist}")
        print(f"  Weights: {weights}")
        
        interpolated_values[i] = np.sum(weights * values_obs[valid_idx]) / np.sum(weights)
    
    return interpolated_values

def write_netcdf(filename, lats, lons, values, variable_name='pm25'):
    with nc.Dataset(filename, 'w', format='NETCDF4') as ds:
        ds.createDimension('lat', lats.shape[0])
        ds.createDimension('lon', lons.shape[1])
        
        lat_var = ds.createVariable('lat', 'f4', ('lat', 'lon'))
        lon_var = ds.createVariable('lon', 'f4', ('lat', 'lon'))
        pm_var = ds.createVariable(variable_name, 'f4', ('lat', 'lon'))
        
        lat_var[:, :] = lats
        lon_var[:, :] = lons
        pm_var[:, :] = values

def plot_surface(lats, lons, values, title, output_file):
    plt.figure(figsize=(12, 8))
    plt.contourf(lons, lats, values, cmap='viridis', levels=np.linspace(np.nanmin(values), np.nanmax(values), 100))
    plt.colorbar(label='PM2.5 ($\mu$g/m$^3$)')
    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(output_file, dpi=300)
    plt.close()

def main():
    # Read input data
    f1 = nc.Dataset('pm25.nc', 'r')
    tlon = f1['MetaData/longitude'][:]
    tlat = f1['MetaData/latitude'][:]
    pm = f1['ObsValue/particulatematter2p5Surface'][:]
    
    # Flatten arrays
    xy_obs = np.column_stack((tlon.flatten(), tlat.flatten()))
    values_obs = pm.flatten()
    
    # Read target grid
    grid_ds = nc.Dataset('fv3_grid_spec.nc', 'r')
    lats = grid_ds.variables['grid_latt'][:]
    lons = grid_ds.variables['grid_lont'][:]
    
    # Flatten grid for interpolation
    xy_grid = np.column_stack((lons.flatten(), lats.flatten()))
    
    # Perform interpolation
    print("Starting interpolation...")
    interpolated_values = inverse_distance_weighting(xy_obs, values_obs, xy_grid, cutoff_radius=0.2, power=2)
    print("Interpolation complete.")
    
    # Reshape the interpolated values back to the 2D grid
    interpolated_values = interpolated_values.reshape(lats.shape)
    
    # Write to NetCDF
    write_netcdf('interpolated_pm25.nc', lats, lons, interpolated_values)
    
    # Plot the surface
    plot_surface(lats, lons, interpolated_values, 'Interpolated PM2.5 Surface', 'interpolated_pm25.png')

if __name__ == '__main__':
    main()


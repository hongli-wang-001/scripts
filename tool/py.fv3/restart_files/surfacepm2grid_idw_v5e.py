import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from matplotlib.colors import Normalize, ListedColormap

def inverse_distance_weighting(xy_obs, values_obs, xy_grid, cutoff_radius, power=2):
    """Perform IDW interpolation with a specified cutoff radius."""
    tree = cKDTree(xy_obs)
    distances, indices = tree.query(xy_grid, k=len(xy_obs), distance_upper_bound=cutoff_radius)
    
    interpolated_values = np.full(len(xy_grid), np.nan)  # Initialize with NaNs
    num_points = len(xy_grid)
    
    for i, (dist, idx) in enumerate(zip(distances, indices)):
        # Progress update
        if i % (num_points // 10) == 0:
            print(f"Interpolating grid point {i}/{num_points}...")
        
        # Filter out infinite distances
        valid_idx = idx[dist < np.inf]
        valid_dist = dist[dist < np.inf]
        
        if valid_idx.size < 3:  # Not enough observations within cutoff radius
            continue
        
        # Calculate weights
        weights = 1 / (valid_dist ** power + np.finfo(float).eps)  # Avoid division by zero
        
        # Debug: Print some intermediate values
        if i < 10:  # Print first 10 for debugging
            print(f"Grid point {i}:")
            print(f"Distances: {valid_dist}")
            print(f"Indices: {valid_idx}")
            print(f"Weights: {weights}")
            print(f"Values: {values_obs[valid_idx]}")
        
        # Interpolate
        interpolated_values[i] = np.sum(weights * values_obs[valid_idx]) / np.sum(weights)
    
    return interpolated_values

def write_netcdf(filename, lats, lons, values, variable_name='pm25'):
    """Write interpolated data to a NetCDF file."""
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
    """Plot the surface and save the figure."""
    plt.figure(figsize=(12, 8))
    
    # Mask NaNs in data
    masked_values = np.ma.masked_invalid(values)
    
    # Determine plot limits
    min_val = np.nanmin(masked_values)
    max_val = np.nanmax(masked_values)
    
    # Check if min and max are valid
    if np.isnan(min_val) or np.isnan(max_val):
        raise ValueError("Cannot determine valid axis limits due to NaNs in data.")
    
    # Plot
    plt.contourf(lons, lats, masked_values, cmap='viridis', levels=np.linspace(min_val, max_val, 100))
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
    
    # Plot original PM2.5 data
    plt.figure(figsize=(12, 8))
    plt.scatter(tlon, tlat, c=pm, cmap='viridis', marker='o', s=10, edgecolor='k')
    plt.colorbar(label='PM2.5 ($\mu$g/m$^3$)')
    plt.title('Original PM2.5 Observations')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig('original_pm25.png', dpi=300)
    plt.close()
    
    # Read target grid
    grid_ds = nc.Dataset('fv3_grid_spec.nc', 'r')
    lats = grid_ds.variables['grid_latt'][:]
    lons = grid_ds.variables['grid_lont'][:]
    
    # Flatten grid for interpolation
    xy_grid = np.column_stack((lons.flatten(), lats.flatten()))
    
    # Perform interpolation with increasing cutoff radius
    max_radius = 1.0
    step = 0.1
    interpolated_values = np.full(lats.shape, np.nan)
    
    for radius in np.arange(0.2, max_radius + step, step):
        print(f"Using cutoff radius: {radius}")
        interpolated_flat = inverse_distance_weighting(xy_obs, values_obs, xy_grid, cutoff_radius=radius, power=2)
        
        # Reshape and check results
        interpolated_values = interpolated_flat.reshape(lats.shape)
        
        # Check if there are valid interpolated values (non-NaN)
        if not np.all(np.isnan(interpolated_values)):
            print("Interpolation complete.")
            break
    
    # Write to NetCDF
    write_netcdf('interpolated_pm25.nc', lats, lons, interpolated_values)
    
    # Plot the surface
    plot_surface(lats, lons, interpolated_values, 'Interpolated PM2.5 Surface', 'interpolated_pm25.png')

if __name__ == '__main__':
    main()


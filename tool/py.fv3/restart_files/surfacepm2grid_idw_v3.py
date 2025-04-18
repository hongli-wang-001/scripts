import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from matplotlib.colors import Normalize, BoundaryNorm
import matplotlib.colors as mcolors

def inverse_distance_weighting(xy_obs, values_obs, xy_grid, cutoff_radius=0.2, power=2):
    """
    Perform Inverse Distance Weighting (IDW) interpolation.
    
    Parameters:
    - xy_obs: (N, 2) array of observed points (latitude, longitude).
    - values_obs: (N,) array of observed values.
    - xy_grid: (M, 2) array of grid points (latitude, longitude).
    - cutoff_radius: Cutoff radius for weighting, beyond which distances are ignored.
    - power: Power parameter for IDW.
    
    Returns:
    - Interpolated values at grid points.
    """
    tree = cKDTree(xy_obs)
    distances, indices = tree.query(xy_grid, k=len(xy_obs), distance_upper_bound=cutoff_radius)
    
    interpolated_values = np.zeros(len(xy_grid))
    for i, (dist, idx) in enumerate(zip(distances, indices)):
        weights = 1 / (dist ** power + np.finfo(float).eps)  # Avoid division by zero
        weights[dist == np.inf] = 0  # Set weights to zero where distance is infinity
        interpolated_values[i] = np.sum(weights * values_obs[idx]) / np.sum(weights)
    
    return interpolated_values

def write_netcdf(filename, lats, lons, values, variable_name='pm25'):
    """
    Write interpolated data to a NetCDF file.
    
    Parameters:
    - filename: Output NetCDF file name.
    - lats: Latitude array.
    - lons: Longitude array.
    - values: Interpolated values.
    - variable_name: Name of the variable to store in the NetCDF file.
    """
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
    """
    Plot the interpolated surface and save the figure.
    
    Parameters:
    - lats: Latitude array.
    - lons: Longitude array.
    - values: Interpolated values.
    - title: Title for the plot.
    - output_file: File name to save the plot.
    """
    plt.figure(figsize=(12, 8))
    plt.contourf(lons, lats, values, cmap='viridis', levels=np.linspace(np.min(values), np.max(values), 100))
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


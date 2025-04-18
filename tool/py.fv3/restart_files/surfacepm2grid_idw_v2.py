import netCDF4 as nc
import numpy as np
from scipy.spatial import cKDTree

# Define the cut-off radius in degrees
CUTOFF_RADIUS = 0.2

def idw_interpolation(lon, lat, values, target_lon, target_lat, cutoff_radius):
    """
    Interpolate values using Inverse Distance Weighting (IDW) with a cut-off radius.
    
    Parameters:
    - lon: Array of longitudes of observations.
    - lat: Array of latitudes of observations.
    - values: Array of observed values.
    - target_lon: Longitude of the target grid point.
    - target_lat: Latitude of the target grid point.
    - cutoff_radius: Cut-off radius for considering observations.
    
    Returns:
    - Interpolated value at the target grid point.
    """
    # Create a KDTree for fast spatial queries
    tree = cKDTree(np.c_[lon, lat])
    
    # Query for points within the cut-off radius
    indices = tree.query_ball_point([target_lon, target_lat], cutoff_radius)
    
    if not indices:
        return np.nan  # Return NaN if no points are within the cut-off radius
    
    # Extract the distances and values of the nearby points
    nearby_lons = lon[indices]
    nearby_lats = lat[indices]
    nearby_values = values[indices]
    
    # Compute distances
    distances = np.sqrt((nearby_lons - target_lon) ** 2 + (nearby_lats - target_lat) ** 2)
    
    # Compute weights: inverse of distance, with a very large weight for zero distance
    weights = 1 / (distances + 1e-10)  # Add a small number to avoid division by zero
    
    # Normalize weights
    weights_sum = np.sum(weights)
    if weights_sum == 0:
        return np.nan  # Return NaN if weights sum to zero (shouldn't happen if there's data)
    
    normalized_weights = weights / weights_sum
    
    # Compute weighted average
    interpolated_value = np.sum(normalized_weights * nearby_values)
    
    return interpolated_value

# Read data from NetCDF file
f1 = nc.Dataset('pm25.nc', 'r')
lon = f1['MetaData/longitude'][:]
lat = f1['MetaData/latitude'][:]
values = f1['ObsValue/particulatematter2p5Surface'][:]

# Define the target grid
grid_ds = nc.Dataset('fv3_grid_spec.nc', 'r')
grid_lats = grid_ds.variables['grid_latt'][:]
grid_lons = grid_ds.variables['grid_lont'][:]

# Create an empty array for interpolated values
interpolated_pm = np.full(grid_lats.shape, np.nan)

# Perform interpolation for each grid point and print progress
total_points = grid_lats.size
current_point = 0

print("Starting interpolation...")

for i in range(grid_lats.shape[0]):
    for j in range(grid_lats.shape[1]):
        interpolated_pm[i, j] = idw_interpolation(lon, lat, values, grid_lons[i, j], grid_lats[i, j], CUTOFF_RADIUS)
        current_point += 1
        
        # Print progress
        if current_point % (total_points // 10) == 0:  # Print progress every 10%
            progress_percentage = (current_point / total_points) * 100
            print(f"Progress: {progress_percentage:.1f}%")

print("Interpolation completed.")

# Save the interpolated data to a NetCDF file or process further


import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import netCDF4 as nc

def great_circle_distance(lat1, lon1, lat2, lon2):
#def great_circle_distance(lon1, lat1, lon2, lat2):
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    # Radius of Earth in kilometers
    R = 6371.0
    return R * c

def find_points_within_radius(query_point, data, radius):
    indices_within_radius = []
    #for i, point in enumerate(data):
    for i in range(data.shape[1]):
        distance = great_circle_distance(query_point[0], query_point[1], data[0,i], data[1,i])
        print(f"find_point: {i} {distance} {query_point[0]} {query_point[1]} {data[0,i]} {data[1,i]}")
        if distance <= radius:
            print(f"found_point: {i} {distance} {query_point[0]} {query_point[1]} {data[0,i]} {data[1,i]}")
            indices_within_radius.append(i)
    return indices_within_radius

def idw_interpolation(lats_obs, lons_obs, values_obs, lats_grid, lons_grid, cutoff_radius_max=500.0, cutoff_step=50.0):
    #num_points = lats_grid.shape[0]
    interpolated_values = np.full(lats_grid.shape, np.nan)
    num_obs_grid = np.zeros(lats_grid.shape)
    final_cutoff_radius = np.full(lats_grid.shape, np.nan)

    # Create the KDTree for observation points
    #tree = cKDTree(np.vstack([lons_obs, lats_obs]).T)

    # Reshape to 1D
    lats_1d = lats_grid.reshape(-1)
    lons_1d = lons_grid.reshape(-1)
    #print("Shape of 1D array:", lats_1d.shape)
    print("Size of 1D array:", lons_1d.size)
    num_points = lats_1d.size

    data = np.array([lats_obs, lons_obs]) 
    print("Shape of 2D array:", data.shape)
    print("Size of 2D array:", data.size)
    print("Size of lats_grid.shape0:", lats_grid.shape[0])
    print("Size of lats_grid.shape1:", lats_grid.shape[1])
    # Iterate over each grid point
    for i in range(lats_grid.shape[0]):  # Loop over rows
        #for j in range(lats_grid.shape[1]):  # Loop over columns
        for j in range(47, 67, 2):
            lat_grid = lats_grid[i, j]
            lon_grid = lons_grid[i, j]
            query_point = np.array([lat_grid, lon_grid])

            # Try different cutoff radii
            for cutoff_radius in np.arange(100.0, cutoff_radius_max + cutoff_step, cutoff_step):
                indices_within_radius = find_points_within_radius(query_point, data, cutoff_radius)
                
                if len(indices_within_radius) > 0:
                    lons_close = lons_obs[indices_within_radius]
                    lats_close = lats_obs[indices_within_radius]
                    values_close = values_obs[indices_within_radius]

                    # Calculate great-circle distances
                    dists = np.array([great_circle_distance(lat_grid, lon_grid, lat, lon) for lat, lon in zip(lats_close, lons_close)])

                    # Calculate weights
                    weights = 1 / np.where(dists == 0, 1e-10, dists)

                    # Compute weighted average
                    weighted_sum = np.sum(weights * values_close)
                    sum_weights = np.sum(weights)

                    if sum_weights > 0:
                        interpolated_values[i, j] = weighted_sum / sum_weights
                        num_obs_grid[i, j] = len(indices_within_radius)
                        final_cutoff_radius[i, j] = cutoff_radius
                        break  # Exit radius loop if we have valid interpolation

    return interpolated_values, num_obs_grid, final_cutoff_radius

def plot_surface(lats, lons, values, title, output_file):
    plt.figure(figsize=(10, 7))
    plt.contourf(lons, lats, values, cmap='viridis', levels=np.linspace(np.nanmin(values), np.nanmax(values), 100))
    plt.colorbar(label='PM2.5 ($\mu$g/m$^3$)')
    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(output_file)
    plt.show()

def write_to_netcdf(lats_grid, lons_grid, interpolated_values, num_obs_grid, final_cutoff_radius, output_file):
    with nc.Dataset(output_file, 'w', format='NETCDF4') as ds:
        # Create dimensions
        ds.createDimension('lat', lats_grid.shape[0])
        ds.createDimension('lon', lons_grid.shape[1])

        # Create variables
        lat_var = ds.createVariable('lat', 'f4', ('lat',))
        lon_var = ds.createVariable('lon', 'f4', ('lon',))
        interp_var = ds.createVariable('interpolated_pm25', 'f4', ('lat', 'lon'))
        num_obs_var = ds.createVariable('num_obs', 'i4', ('lat', 'lon'))
        cutoff_var = ds.createVariable('cutoff_radius', 'f4', ('lat', 'lon'))

        # Assign data to variables
        lat_var[:] = lats_grid[:, 0]  # Assuming lats_grid is 2D with lat values in the first column
        lon_var[:] = lons_grid[0, :]  # Assuming lons_grid is 2D with lon values in the first row
        interp_var[:, :] = interpolated_values
        num_obs_var[:, :] = num_obs_grid
        cutoff_var[:, :] = final_cutoff_radius

        # Add attributes
        lat_var.units = 'degrees_north'
        lon_var.units = 'degrees_east'
        interp_var.units = 'µg/m³'
        num_obs_var.units = 'count'
        cutoff_var.units = 'degrees'
        interp_var.description = 'Interpolated PM2.5 values'
        num_obs_var.description = 'Number of observations used for interpolation'
        cutoff_var.description = 'Final cutoff radius used for interpolation'

def main():
    # Load observation data
    f1 = nc.Dataset('pm25.nc', 'r')
    lons_obs = f1['MetaData/longitude'][:]
    lats_obs = f1['MetaData/latitude'][:]
    values_obs = f1['ObsValue/particulatematter2p5Surface'][:]
    
    # Load grid data
    grid_ds = nc.Dataset('fv3_grid_spec.nc', 'r')
    lats_grid = grid_ds.variables['grid_latt'][:]
    lons_grid = grid_ds.variables['grid_lont'][:]

    # Perform interpolation
    interpolated_values, num_obs_grid, final_cutoff_radius = idw_interpolation(
        lats_obs, lons_obs, values_obs, lats_grid, lons_grid
    )
    
    # Plot the results
    plot_surface(lats_grid, lons_grid, interpolated_values, 'Interpolated PM2.5 Surface', 'interpolated_pm25_x47.png')

    # Write results to NetCDF
    write_to_netcdf(lats_grid, lons_grid, interpolated_values, num_obs_grid, final_cutoff_radius, 'interpolated_pm25_x47.nc')

if __name__ == "__main__":
    main()


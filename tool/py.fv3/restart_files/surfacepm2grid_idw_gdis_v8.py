import numpy as np
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import netCDF4 as nc

def great_circle_distance(lat1, lon1, lat2, lon2):
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    R = 6371.0
    return R * c

def idw_interpolation(lats_obs, lons_obs, values_obs, lats_grid, lons_grid, cutoff_radius_max=15.0, cutoff_step=10.0):
    interpolated_values = np.full(lats_grid.shape, np.nan)
    num_obs_grid = np.zeros(lats_grid.shape)
    final_cutoff_radius = np.full(lats_grid.shape, np.nan)

    # Create the KDTree for observation points
    tree = cKDTree(np.vstack([lons_obs, lats_obs]).T)

    # Reshape grid for easy looping
    lats_1d = lats_grid.flatten()
    lons_1d = lons_grid.flatten()
    
    # Iterate over each grid point
    for idx in range(len(lats_1d)):
        lat_grid = lats_1d[idx]
        lon_grid = lons_1d[idx]
        query_point = np.array([[lon_grid, lat_grid]])

        # Find indices within varying radii using cKDTree
        for cutoff_radius in np.arange(10.0, cutoff_radius_max + cutoff_step, cutoff_step):
            indices_within_radius = tree.query_ball_point(query_point, cutoff_radius)

            if indices_within_radius:
                lons_close = lons_obs[indices_within_radius]
                lats_close = lats_obs[indices_within_radius]
                values_close = values_obs[indices_within_radius]

                # Calculate distances
                dists = great_circle_distance(lat_grid, lon_grid, lats_close, lons_close)
                weights = 1 / np.where(dists == 0, 1e-10, dists)

                # Compute weighted average
                weighted_sum = np.sum(weights * values_close)
                sum_weights = np.sum(weights)

                if sum_weights > 0:
                    interpolated_values[idx // lats_grid.shape[1], idx % lats_grid.shape[1]] = weighted_sum / sum_weights
                    num_obs_grid[idx // lats_grid.shape[1], idx % lats_grid.shape[1]] = len(indices_within_radius)
                    final_cutoff_radius[idx // lats_grid.shape[1], idx % lats_grid.shape[1]] = cutoff_radius
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
        ds.createDimension('lat', lats_grid.shape[0])
        ds.createDimension('lon', lons_grid.shape[1])

        lat_var = ds.createVariable('lat', 'f4', ('lat',))
        lon_var = ds.createVariable('lon', 'f4', ('lon',))
        interp_var = ds.createVariable('interpolated_pm25', 'f4', ('lat', 'lon'))
        num_obs_var = ds.createVariable('num_obs', 'i4', ('lat', 'lon'))
        cutoff_var = ds.createVariable('cutoff_radius', 'f4', ('lat', 'lon'))

        lat_var[:] = lats_grid[:, 0]
        lon_var[:] = lons_grid[0, :]
        interp_var[:, :] = interpolated_values
        num_obs_var[:, :] = num_obs_grid
        cutoff_var[:, :] = final_cutoff_radius

        lat_var.units = 'degrees_north'
        lon_var.units = 'degrees_east'
        interp_var.units = 'µg/m³'
        num_obs_var.units = 'count'
        cutoff_var.units = 'degrees'
        interp_var.description = 'Interpolated PM2.5 values'
        num_obs_var.description = 'Number of observations used for interpolation'
        cutoff_var.description = 'Final cutoff radius used for interpolation'

def main():
    f1 = nc.Dataset('pm25.nc', 'r')
    lons_obs = f1['MetaData/longitude'][:]
    lats_obs = f1['MetaData/latitude'][:]
    values_obs = f1['ObsValue/particulatematter2p5Surface'][:]

    grid_ds = nc.Dataset('fv3_grid_spec.nc', 'r')
    lats_grid = grid_ds.variables['grid_latt'][:]
    lons_grid = grid_ds.variables['grid_lont'][:]

    interpolated_values, num_obs_grid, final_cutoff_radius = idw_interpolation(
        lats_obs, lons_obs, values_obs, lats_grid, lons_grid
    )

    plot_surface(lats_grid, lons_grid, interpolated_values, 'Interpolated PM2.5 Surface', 'interpolated_pm25.png')
    write_to_netcdf(lats_grid, lons_grid, interpolated_values, num_obs_grid, final_cutoff_radius, 'interpolated_pm25.nc')

if __name__ == "__main__":
    main()


import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree

def idw_interpolation(lats_obs, lons_obs, values_obs, lats_grid, lons_grid, cutoff_radius_max=1.0, cutoff_step=0.1):
    num_points = lats_grid.shape[0]
    interpolated_values = np.full(lats_grid.shape, np.nan)
    num_obs_grid = np.zeros(lats_grid.shape)
    final_cutoff_radius = np.full(lats_grid.shape, np.nan)

    tree = cKDTree(np.vstack([lons_obs, lats_obs]).T)

    for i in range(num_points):
        lat_grid = lats_grid[i]
        lon_grid = lons_grid[i]
        point = np.array([lon_grid, lat_grid])  # Ensure this is a 1D array

        for cutoff_radius in np.arange(0.2, cutoff_radius_max + cutoff_step, cutoff_step):
            idx = tree.query_ball_point(point, cutoff_radius)
            if len(idx) > 0:
                lons_close = lons_obs[idx]
                lats_close = lats_obs[idx]
                values_close = values_obs[idx]

                dists = np.sqrt((lons_close - lon_grid)**2 + (lats_close - lat_grid)**2)
                weights = 1 / np.where(dists == 0, 1e-10, dists)

                weighted_sum = np.sum(weights * values_close)
                sum_weights = np.sum(weights)

                if sum_weights > 0:
                    interpolated_values[i] = weighted_sum / sum_weights
                    num_obs_grid[i] = len(idx)
                    final_cutoff_radius[i] = cutoff_radius
                    break

    return interpolated_values, num_obs_grid, final_cutoff_radius

def save_to_netcdf(filename, lats_grid, lons_grid, interpolated_values, num_obs_grid, final_cutoff_radius):
    with nc.Dataset(filename, 'w', format='NETCDF4') as ds:
        ds.createDimension('lat', lats_grid.shape[0])
        ds.createDimension('lon', lons_grid.shape[1])
        
        latitudes = ds.createVariable('latitude', 'f4', ('lat',))
        longitudes = ds.createVariable('longitude', 'f4', ('lon',))
        pm2_5 = ds.createVariable('pm25', 'f4', ('lat', 'lon',))
        num_obs = ds.createVariable('num_observations', 'i4', ('lat', 'lon',))
        cutoff_radius = ds.createVariable('cutoff_radius', 'f4', ('lat', 'lon',))
        
        latitudes[:] = lats_grid[:, 0]  # Assuming latitudes are in the first column
        longitudes[:] = lons_grid[0, :]  # Assuming longitudes are in the first row
        pm2_5[:, :] = interpolated_values
        num_obs[:, :] = num_obs_grid
        cutoff_radius[:, :] = final_cutoff_radius

def plot_surface(lats, lons, values, title, filename):
    plt.figure(figsize=(12, 6))
    plt.contourf(lons, lats, values, cmap='viridis', levels=np.linspace(np.nanmin(values), np.nanmax(values), 100))
    plt.colorbar(label='PM2.5 ($\mu$g/m$^3$)')
    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(filename, dpi=300)
    plt.close()

def main():
    # Load observational data
    obs_file = 'pm25.nc'
    with nc.Dataset(obs_file, 'r') as f:
        lats_obs = f['MetaData/latitude'][:]
        lons_obs = f['MetaData/longitude'][:]
        values_obs = f['ObsValue/particulatematter2p5Surface'][:]
    
    # Load grid data
    grid_file = 'fv3_grid_spec.nc'
    with nc.Dataset(grid_file, 'r') as grid_ds:
        lats_grid = grid_ds.variables['grid_latt'][:]
        lons_grid = grid_ds.variables['grid_lont'][:]
    
    # Perform IDW interpolation
    interpolated_values, num_obs_grid, final_cutoff_radius = idw_interpolation(
        lats_obs, lons_obs, values_obs, lats_grid, lons_grid, cutoff_radius_max=1.0, cutoff_step=0.1
    )
    
    # Save results to NetCDF
    output_nc_file = 'interpolated_pm25.nc'
    save_to_netcdf(output_nc_file, lats_grid, lons_grid, interpolated_values, num_obs_grid, final_cutoff_radius)
    
    # Plot the interpolated surface
    plot_surface(lats_grid, lons_grid, interpolated_values, 'Interpolated PM2.5 Surface', 'interpolated_pm25.png')
    
    print("Interpolation and plotting complete.")

if __name__ == "__main__":
    main()


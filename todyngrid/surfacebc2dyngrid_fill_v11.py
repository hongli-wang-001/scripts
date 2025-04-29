import numpy as np
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

def update_grid_with_observations(lats_obs, lons_obs, values_obs, lats_grid, lons_grid, grid_size=0.1):
    # Convert lons_obs from [-180, 180] to [0, 360]
    lons_obs = np.where(lons_obs < 0, lons_obs + 360, lons_obs)

    interpolated_values = np.full(lats_grid.shape, np.nan)

    for obs_idx in range(len(lats_obs)):
        lat_obs = lats_obs[obs_idx]
        lon_obs = lons_obs[obs_idx]
        value_obs = values_obs[obs_idx]

        lat_min = lat_obs - 15 / 111.32
        lat_max = lat_obs + 15 / 111.32
        lon_min = lon_obs - 15 / (111.32 * np.cos(np.radians(lat_obs)))
        lon_max = lon_obs + 15 / (111.32 * np.cos(np.radians(lat_obs)))

        print(f'Observation {obs_idx}, {lat_obs}, {lon_obs}: Lat Range: ({lat_min}, {lat_max}), Lon Range: ({lon_min}, {lon_max})')

        #lat_indices = np.where((lats_grid >= lat_min) & (lats_grid <= lat_max))[0]
        #lon_indices = np.where((lons_grid >= lon_min) & (lons_grid <= lon_max))[1]
        # Create a 2D mask for both latitude and longitude
        lat_mask = (lats_grid >= lat_min) & (lats_grid <= lat_max)
        lon_mask = (lons_grid >= lon_min) & (lons_grid <= lon_max)
        combined_mask = lat_mask & lon_mask

        # Get indices where both masks are true
        lat_indices, lon_indices = np.where(combined_mask)

        found_valid_point = False

        #for i in lat_indices:
        #    for j in lon_indices:
        for i, j in zip(lat_indices, lon_indices):
                distance = great_circle_distance(lat_obs, lon_obs, lats_grid[i, j], lons_grid[i, j])
                print(f'Distance to Grid Point (Lat: {lats_grid[i, j]}, Lon: {lons_grid[i, j]}): {distance:.2f} km')
                print(f'Observation at (Lat: {lat_obs}, Lon: {lon_obs}, Value: {value_obs}) updates Grid Point at (Lat: {lats_grid[i, j]}, Lon: {lons_grid[i, j]})')
                if np.isnan(interpolated_values[i, j]):
                    interpolated_values[i, j] = value_obs
                else:
                    interpolated_values[i, j] = max(interpolated_values[i, j], value_obs)
                    found_valid_point = True

        #if not found_valid_point:
        #    print(f'No valid grid point found for observation {obs_idx} at (Lat: {lat_obs}, Lon: {lon_obs})')

        if obs_idx % 100 == 0:
            print(f'Processed {obs_idx} observations out of {len(lats_obs)}')

    return interpolated_values

def update_grid_with_observations_2(lats_obs, lons_obs, values_obs, lats_grid, lons_grid, grid_size=0.1):
    # Convert lons_obs from [-180, 180] to [0, 360]
    lons_obs = np.where(lons_obs < 0, lons_obs + 360, lons_obs)

    interpolated_values = np.full(lats_grid.shape, np.nan)

    for obs_idx in range(len(lats_obs)):
        lat_obs = lats_obs[obs_idx]
        lon_obs = lons_obs[obs_idx]
        value_obs = values_obs[obs_idx]
        if np.isnan(value_obs):
           continue  # Skip this observation if value is NaN
        if value_obs == -1:
           continue  # Skip this observation if value is -1 (missing)


        lat_min = lat_obs - 15 / 111.32
        lat_max = lat_obs + 15 / 111.32
        lon_min = lon_obs - 15 / (111.32 * np.cos(np.radians(lat_obs)))
        lon_max = lon_obs + 15 / (111.32 * np.cos(np.radians(lat_obs)))

        print(f'Observation {obs_idx}, {lat_obs}, {lon_obs}: Lat Range: ({lat_min}, {lat_max}), Lon Range: ({lon_min}, {lon_max})')

        #lat_indices = np.where((lats_grid >= lat_min) & (lats_grid <= lat_max))[0]
        #lon_indices = np.where((lons_grid >= lon_min) & (lons_grid <= lon_max))[1]
        # Create a 2D mask for both latitude and longitude
        lat_mask = (lats_grid >= lat_min) & (lats_grid <= lat_max)
        lon_mask = (lons_grid >= lon_min) & (lons_grid <= lon_max)
        combined_mask = lat_mask & lon_mask

        # Get indices where both masks are true
        lat_indices, lon_indices = np.where(combined_mask)

        found_valid_point = False
        dis_temp = 100
        #for i in lat_indices:
        #    for j in lon_indices:
        for i, j in zip(lat_indices, lon_indices):
                distance = great_circle_distance(lat_obs, lon_obs, lats_grid[i, j], lons_grid[i, j])
                print(f'Distance to Grid Point (Lat: {lats_grid[i, j]}, Lon: {lons_grid[i, j]}): {distance:.2f} km')
                if distance <= dis_temp:
                   dis_temp = distance
                   ni , nj = i , j
                   print(f'Observation at (Lat: {lat_obs}, Lon: {lon_obs}, Dis: {dis_temp}')
        if dis_temp <= 13:
                    print(f'Observation at (Lat: {lat_obs}, Lon: {lon_obs}, Value: {value_obs}) updates Grid Point at (Lat: {lats_grid[i, j]}, Lon: {lons_grid[i, j]})')
                    
                    #for di in [-1, 0, 1]:
                    #    for dj in [-1, 0, 1]:
                    #        ni, nj = i + di, j + dj
                    if 0 <= ni < lats_grid.shape[0] and 0 <= nj < lons_grid.shape[1]:
                                if np.isnan(interpolated_values[ni, nj]):
                                    interpolated_values[ni, nj] = value_obs
                                #else:
                                #    interpolated_values[i, j] = max(interpolated_values[i, j], value_obs)
                    found_valid_point = True
                    #break

        if not found_valid_point:
            print(f'No valid grid point found for observation {obs_idx} at (Lat: {lat_obs}, Lon: {lon_obs})')

        if obs_idx % 100 == 0:
            print(f'Processed {obs_idx} observations out of {len(lats_obs)}')

    return interpolated_values

def update_grid_with_observations_3(lats_obs, lons_obs, values_obs, lats_grid, lons_grid, grid_size=0.1):
    # Convert lons_obs from [-180, 180] to [0, 360]
    lons_obs = np.where(lons_obs < 0, lons_obs + 360, lons_obs)

    interpolated_values = np.full(lats_grid.shape, np.nan)

    for obs_idx in range(len(lats_obs)):
        lat_obs = lats_obs[obs_idx]
        lon_obs = lons_obs[obs_idx]
        value_obs = values_obs[obs_idx]

        lat_min = lat_obs - 15 / 111.32
        lat_max = lat_obs + 15 / 111.32
        lon_min = lon_obs - 15 / (111.32 * np.cos(np.radians(lat_obs)))
        lon_max = lon_obs + 15 / (111.32 * np.cos(np.radians(lat_obs)))

        print(f'Observation {obs_idx}, {lat_obs}, {lon_obs}: Lat Range: ({lat_min}, {lat_max}), Lon Range: ({lon_min}, {lon_max})')

        #lat_indices = np.where((lats_grid >= lat_min) & (lats_grid <= lat_max))[0]
        #lon_indices = np.where((lons_grid >= lon_min) & (lons_grid <= lon_max))[1]
        # Create a 2D mask for both latitude and longitude
        lat_mask = (lats_grid >= lat_min) & (lats_grid <= lat_max)
        lon_mask = (lons_grid >= lon_min) & (lons_grid <= lon_max)
        combined_mask = lat_mask & lon_mask

        # Get indices where both masks are true
        lat_indices, lon_indices = np.where(combined_mask)

        found_valid_point = False

        #for i in lat_indices:
        #    for j in lon_indices:
        for i, j in zip(lat_indices, lon_indices):
                distance = great_circle_distance(lat_obs, lon_obs, lats_grid[i, j], lons_grid[i, j])
                print(f'Distance to Grid Point (Lat: {lats_grid[i, j]}, Lon: {lons_grid[i, j]}): {distance:.2f} km')

                if distance <= 15:
                    print(f'Observation at (Lat: {lat_obs}, Lon: {lon_obs}, Value: {value_obs}) updates Grid Point at (Lat: {lats_grid[i, j]}, Lon: {lons_grid[i, j]})')

                    for di in [-1, 0, 1]:
                        for dj in [-1, 0, 1]:
                            ni, nj = i + di, j + dj
                            if 0 <= ni < lats_grid.shape[0] and 0 <= nj < lons_grid.shape[1]:
                                if np.isnan(interpolated_values[ni, nj]):
                                    interpolated_values[ni, nj] = value_obs
                                else:
                                    interpolated_values[ni, nj] = max(interpolated_values[ni, nj], value_obs)
                    found_valid_point = True
                    break

        if not found_valid_point:
            print(f'No valid grid point found for observation {obs_idx} at (Lat: {lat_obs}, Lon: {lon_obs})')

        if obs_idx % 100 == 0:
            print(f'Processed {obs_idx} observations out of {len(lats_obs)}')

    return interpolated_values

def plot_surface(lats, lons, values, title, output_file):
    if np.all(np.isnan(values)):
        print("No valid data available for plotting.")
        return

    values_filled = np.nan_to_num(values, nan=0)

    plt.figure(figsize=(10, 7))
    plt.contourf(lons, lats, values_filled, cmap='viridis', levels=np.linspace(np.nanmin(values_filled), np.nanmax(values_filled), 100))
    plt.colorbar(label='PM2.5 ($\mu$g/m$^3$)')
    plt.title(title)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.savefig(output_file)
    plt.show()

def write_to_netcdf_ll(lats_grid, lons_grid, interpolated_values, output_file):
    with nc.Dataset(output_file, 'w', format='NETCDF4') as ds:
        ds.createDimension('lat', lats_grid.shape[0])
        ds.createDimension('lon', lons_grid.shape[1])
        lat_var = ds.createVariable('lat', 'f4', ('lat',))
        lon_var = ds.createVariable('lon', 'f4', ('lon',))
        interp_var = ds.createVariable('interpolated_pm25', 'f4', ('lat', 'lon'))

        lat_var[:] = lats_grid[:, 0]
        lon_var[:] = lons_grid[0, :]
        interp_var[:, :] = interpolated_values

        lat_var.units = 'degrees_north'
        lon_var.units = 'degrees_east'
        interp_var.units = 'µg/m³'
        interp_var.description = 'Interpolated PM2.5 values'

def write_to_netcdf(lats_grid, lons_grid, interpolated_values, output_file):
    with nc.Dataset(output_file, 'w', format='NETCDF4') as ds:
        # Dimensions match the 2D grid
        ny, nx = lats_grid.shape
        ds.createDimension('y', ny)
        ds.createDimension('x', nx)

        # Create 2D lat/lon variables
        lat_var = ds.createVariable('latitude', 'f4', ('y', 'x'))
        lon_var = ds.createVariable('longitude', 'f4', ('y', 'x'))

        # Create interpolated data variable
        interp_var = ds.createVariable('bc', 'f4', ('y', 'x'), zlib=True, complevel=4)

        # Assign data to variables
        lat_var[:, :] = lats_grid
        lon_var[:, :] = lons_grid
        interp_var[:, :] = interpolated_values

        # Add variable metadata
        lat_var.units = 'degrees_north'
        lon_var.units = 'degrees_east'
        interp_var.units = 'µg/m³'
        interp_var.long_name = 'BC concentration'

        # Optionally: Add global attributes
        ds.title = 'PM2.5 Interpolation Output'
        ds.projection = 'Lambert Conformal Conic'

def main():
    #f1 = nc.Dataset('pm25.nc', 'r')
    #lons_obs = f1['MetaData/longitude'][:]
    #lats_obs = f1['MetaData/latitude'][:]
    #values_obs = f1['ObsValue/particulatematter2p5Surface'][:]
    f1 = nc.Dataset('pm25.nc', 'r')
    lons_obs = f1['longitude'][:].flatten()
    lats_obs = f1['latitude'][:].flatten()
    values_obs = f1['BC'][:].flatten()  # assuming 'BC' is 2D or 1xN

    grid_ds = nc.Dataset('dyn_grid_spec.nc', 'r')
    lats_grid = grid_ds.variables['lat'][:]
    lons_grid = grid_ds.variables['lon'][:]

    # Print statistics for obs values
    print("Obs Values:")
    print(f"Min: {np.nanmin(values_obs)}, Max: {np.nanmax(values_obs)}, Mean: {np.nanmean(values_obs)}, Std: {np.nanstd(values_obs)}")
    interpolated_values = update_grid_with_observations_2(lats_obs, lons_obs, values_obs, lats_grid, lons_grid)
    # Print statistics for interpolated values
    print("Interpolated Values:")
    print(f"Min: {np.nanmin(interpolated_values)}, Max: {np.nanmax(interpolated_values)}, Mean: {np.nanmean(interpolated_values)}, Std: {np.nanstd(interpolated_values)}")

    num_valid_obs = np.count_nonzero(~np.isnan(values_obs) & (values_obs != -1))
    num_valid_points = np.count_nonzero(~np.isnan(interpolated_values))
    print(f'Number of valid interpolated points: {num_valid_points}, valid_obs: {num_valid_obs}')

    plot_surface(lats_grid, lons_grid, interpolated_values, 'Interpolated PM2.5 Surface', 'v10_interpolated_pm25.png')

    write_to_netcdf(lats_grid, lons_grid, interpolated_values, 'v10_interpolated_pm25.nc')

if __name__ == "__main__":
    main()


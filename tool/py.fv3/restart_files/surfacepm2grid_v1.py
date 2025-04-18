#!/usr/bin/env python
import netCDF4 as nc
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

# Load the surface PM2.5 data
input_file = 'pm25.nc'
f1 = nc.Dataset(input_file, 'r')
tlon = f1['MetaData/longitude'][:]
tlat = f1['MetaData/latitude'][:]
pm = f1['ObsValue/particulatematter2p5Surface'][:].flatten()  # Flatten if necessary

# Load the target grid (2D irregular grid)
grid_file = 'fv3_grid_spec.nc'
grid_ds = nc.Dataset(grid_file, 'r')
lats = grid_ds.variables['grid_latt'][:]
lons = grid_ds.variables['grid_lont'][:]

# Flatten the target grid to make it compatible with griddata
lats_flat = lats.flatten()
lons_flat = lons.flatten()

# Define the target grid
lon_target, lat_target = np.meshgrid(np.unique(lons), np.unique(lats))

# Interpolation using griddata
# Input data for interpolation
points = np.vstack((tlon, tlat)).T  # (N, 2) array of points
values = pm  # (N,) array of PM2.5 values

# Target grid points for interpolation
grid_points = np.vstack((lon_target.flatten(), lat_target.flatten())).T

# Perform the interpolation
pm_interpolated = griddata(points, values, grid_points, method='linear')

# Reshape the interpolated data back to the target grid shape
pm_interpolated_grid = pm_interpolated.reshape(lats.shape)

# Plotting the interpolated data
plt.figure(figsize=(12, 8))
plt.contourf(lons, lats, pm_interpolated_grid, cmap='viridis', levels=100)
plt.colorbar(label='PM2.5 Concentration')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Interpolated PM2.5 Data on Irregular Grid')
plt.show()

# Optionally, save the interpolated data to a new NetCDF file
output_file = 'interpolated_pm25.nc'
with nc.Dataset(output_file, 'w', format='NETCDF4') as nc_out:
    nc_out.createDimension('lat', lats.shape[0])
    nc_out.createDimension('lon', lons.shape[1])
    nc_out.createVariable('lat', 'f4', ('lat', 'lon'))
    nc_out.createVariable('lon', 'f4', ('lat', 'lon'))
    nc_out.createVariable('pm25', 'f4', ('lat', 'lon'))

    nc_out.variables['lat'][:] = lats
    nc_out.variables['lon'][:] = lons
    nc_out.variables['pm25'][:] = pm_interpolated_grid


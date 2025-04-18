import sys
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy.interpolate import griddata
from cartopy import crs as ccrs
import netCDF4 as nc
from pyproj import Proj
import cartopy.feature as cfeature

# Define Lambert conformal projection parameters
lat1 = 22.57418  # Standard parallel 1
lat2 = 38.5      # Standard parallel 2
cen_lat = 38.5
cen_lon = -97.5

# Load CMAQBC data in phy netcdf
bc_file = 'file1.nc'
with nc.Dataset(bc_file, 'r') as bc_nc:
    bc_values = bc_nc.variables['aecj'][0, 63, :, :]
    aeci_values = bc_nc.variables['aeci'][0, 63, :, :]
    bc_values = bc_values + aeci_values

# Load CMAQ lat lon 
ll_file = 'lonlat.nc'
with nc.Dataset(ll_file, 'r') as ll_nc:
    lon = ll_nc.variables['grid_lont'][:]
    lat = ll_nc.variables['grid_latt'][:]
#print(f"lon shape: {lon.shape}, size: {lon.size}")

# Load dry_air_density
dens_file = 'dens.nc'
with nc.Dataset(dens_file, 'r') as dens_nc:
    air_dens = dens_nc.variables['dry_air_density'][0, 63, :, :]

# CMAQ BC ug/Kg to ug/M^3
bc_values = bc_values*air_dens

# Load AIRNOW data
airnow_ds = xr.open_dataset('file2.nc')
obs_bc = airnow_ds['BC'][0, 0].values
latitude = airnow_ds['latitude'][0].values
longitude = airnow_ds['longitude'][0].values
obs_bc[obs_bc == -1] = np.nan  # Replace missing values with NaN

longitude = np.where(longitude < 0, longitude + 360, longitude)

print("min and max (longitude):")
#print(longitude)
print(longitude.min())
print(longitude.max())

# Print shape and size of obs_bc
print(f"obs_bc shape: {obs_bc.shape}, size: {obs_bc.size}")
print(f"latitude shape: {latitude.shape}, size: {latitude.size}")
#print("obs_bc:")
#print(bc_values)

#sys.exit("Exiting the script.")

# Create an array mask 
#nan_mask = np.isnan(obs_bc)
valid_obs_mask = np.isfinite(obs_bc)
# Count the number of True values in valid_mask
num_valid_obs = np.sum(valid_obs_mask)

# Print the number of valid observations
print(f"Number of valid observations: {num_valid_obs}")

# Filter valid observations
valid_longitude = longitude[valid_obs_mask]
valid_latitude = latitude[valid_obs_mask]
valid_obs_bc = obs_bc[valid_obs_mask]


# Interpolation to observation locations
interpolated_bc_a = np.full_like(obs_bc, np.nan)
if longitude.size > 0:
    interpolated_bc_a = griddata(
        (lon.flatten(), lat.flatten()),
        bc_values.flatten(),
        (longitude, latitude),
        method='linear',
        fill_value=np.nan
    )

# Interpolation to observation locations
interpolated_bc = np.full_like(valid_obs_bc, np.nan)
if valid_obs_bc.size > 0:
    interpolated_bc = griddata(
        (lon.flatten(), lat.flatten()),
        bc_values.flatten(),
        (valid_longitude, valid_latitude),
        method='linear',
        fill_value=np.nan
    )
#print("Valid interpolated_bc:")
#print(interpolated_bc)

# Output statistics
print("CMAQ BC Statistics:")
print(f"Mean: {np.nanmean(bc_values):.3f}, Std: {np.nanstd(bc_values):.3f}, "
      f"Min: {np.nanmin(bc_values):.3f}, Max: {np.nanmax(bc_values):.3f}")

print("AIRNOW BC Statistics:")
print(f"Mean: {np.nanmean(valid_obs_bc):.3f}, Std: {np.nanstd(valid_obs_bc):.3f}, "
      f"Min: {np.nanmin(valid_obs_bc):.3f}, Max: {np.nanmax(valid_obs_bc):.3f}")

print("Interpolated BC Statistics:")
print(f"Mean: {np.nanmean(interpolated_bc):.3f}, Std: {np.nanstd(interpolated_bc):.3f}, "
      f"Min: {np.nanmin(interpolated_bc):.3f}, Max: {np.nanmax(interpolated_bc):.3f}")

# Calculate and print differences
bc_difference = interpolated_bc - valid_obs_bc
print("Difference Statistics:")
print(f"Mean: {np.nanmean(bc_difference):.3f}, Std: {np.nanstd(bc_difference):.3f}, "
      f"Min: {np.nanmin(bc_difference):.3f}, Max: {np.nanmax(bc_difference):.3f}")

with nc.Dataset('bc_airnowexcnon_cmaq_interpolated.nc', 'w', format='NETCDF4') as nc_out:
    nc_out.createDimension('points', len(valid_obs_bc))
    lon_out = nc_out.createVariable('longitude', np.float32, ('points',))
    lon_out[:] = valid_longitude
    lat_out = nc_out.createVariable('latitude', np.float32, ('points',))
    lat_out[:] = valid_latitude
    bc_out = nc_out.createVariable('interpolated_bc', np.float32, ('points',))
    bc_out[:] = interpolated_bc
    obs_out = nc_out.createVariable('airnow_bc', np.float32, ('points',))
    obs_out[:] = valid_obs_bc

# Create a Lambert Conformal projection for plotting
lambert_proj = ccrs.LambertConformal(
    central_latitude=cen_lat,
    central_longitude=cen_lon,
    standard_parallels=(lat1, lat2)
)

# Plot the interpolated BC
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

# Initialize arrays for transformed coordinates
longitude_lcc = np.zeros(valid_longitude.shape)
latitude_lcc = np.zeros(valid_latitude.shape)

# Transform the data point to projection coordinates
# Transform each point
for i in range(len(valid_longitude)):
    longitude_lcc[i], latitude_lcc[i] = lambert_proj.transform_point(valid_longitude[i], valid_latitude[i], ccrs.PlateCarree())

scatter = ax.scatter(longitude_lcc, latitude_lcc, c=interpolated_bc, cmap='YlOrRd', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='Interpolated BC')
plt.title('BC Interpolated Values')
plt.savefig('bc_cmaq_interpolated_excnan.png')
plt.show()

# Plot the AIRNOW BC
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

scatter = ax.scatter(longitude_lcc, latitude_lcc, c=valid_obs_bc, cmap='YlOrRd', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='AIRNOW BC')
plt.title('BC AIRNOW Values')
plt.savefig('bc_airnow_excnan.png')
plt.show()

# Plot the difference
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

scatter = ax.scatter(longitude_lcc, latitude_lcc, c=bc_difference, cmap='RdBu_r', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='Difference (Interpolated BC - AIRNOW BC)')
plt.title('Difference between Interpolated BC and AIRNOW BC')
plt.savefig('bc_difference_excnan.png')
plt.show()


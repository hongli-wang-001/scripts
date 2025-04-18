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

# Load CMAQAOD data in phy netcdf
aod_file = 'file1.nc'
with nc.Dataset(aod_file, 'r') as aod_nc:
    aod_values = aod_nc.variables['aod'][0, :, :]

# Load CMAQ lat lon in dyn file
ll_file = 'lonlat.nc'
with nc.Dataset(ll_file, 'r') as ll_nc:
    lon = ll_nc.variables['lon'][:]
    lat = ll_nc.variables['lat'][:]
#print(f"lon shape: {lon.shape}, size: {lon.size}")

# Load AERONET data
aeronet_ds = xr.open_dataset('file2.nc')
obs_aod = aeronet_ds['aod_550nm'][0, 0].values
latitude = aeronet_ds['latitude'].values
longitude = aeronet_ds['longitude'].values
obs_aod[obs_aod == -1] = np.nan  # Replace missing values with NaN

longitude = np.where(longitude < 0, longitude + 360, longitude)

print("min and max (longitude):")
#print(longitude)
print(longitude.min())
print(longitude.max())

# Print shape and size of obs_aod
print(f"obs_aod shape: {obs_aod.shape}, size: {obs_aod.size}")
print(f"latitude shape: {latitude.shape}, size: {latitude.size}")
#print("obs_aod:")
#print(aod_values)

#sys.exit("Exiting the script.")

# Create an array mask 
#nan_mask = np.isnan(obs_aod)
valid_obs_mask = np.isfinite(obs_aod)
# Count the number of True values in valid_mask
num_valid_obs = np.sum(valid_obs_mask)

# Print the number of valid observations
print(f"Number of valid observations: {num_valid_obs}")

# Filter valid observations
valid_longitude = longitude[valid_obs_mask]
valid_latitude = latitude[valid_obs_mask]
valid_obs_aod = obs_aod[valid_obs_mask]


# Interpolation to observation locations
interpolated_aod = np.full_like(obs_aod, np.nan)
if longitude.size > 0:
    interpolated_aod = griddata(
        (lon.flatten(), lat.flatten()),
        aod_values.flatten(),
        (longitude, latitude),
        method='linear',
        fill_value=np.nan
    )

# Interpolation to observation locations
interpolated_aod_b = np.full_like(valid_obs_aod, np.nan)
if valid_longitude.size > 0:
    interpolated_aod_b = griddata(
        (lon.flatten(), lat.flatten()),
        aod_values.flatten(),
        (valid_longitude, valid_latitude),
        method='linear',
        fill_value=np.nan
    )
#print("Valid interpolated_aod:")
#print(interpolated_aod)

# Output statistics
print("CMAQ AOD Statistics:")
print(f"Mean: {np.nanmean(aod_values):.3f}, Std: {np.nanstd(aod_values):.3f}, "
      f"Min: {np.nanmin(aod_values):.3f}, Max: {np.nanmax(aod_values):.3f}")

print("AERONET AOD Statistics:")
print(f"Mean: {np.nanmean(obs_aod):.3f}, Std: {np.nanstd(obs_aod):.3f}, "
      f"Min: {np.nanmin(obs_aod):.3f}, Max: {np.nanmax(obs_aod):.3f}")

print("Interpolated AOD Statistics:")
print(f"Mean: {np.nanmean(interpolated_aod):.3f}, Std: {np.nanstd(interpolated_aod):.3f}, "
      f"Min: {np.nanmin(interpolated_aod):.3f}, Max: {np.nanmax(interpolated_aod):.3f}")

# Calculate and print differences
aod_difference = interpolated_aod - obs_aod
print("Difference Statistics:")
print(f"Mean: {np.nanmean(aod_difference):.3f}, Std: {np.nanstd(aod_difference):.3f}, "
      f"Min: {np.nanmin(aod_difference):.3f}, Max: {np.nanmax(aod_difference):.3f}")

with nc.Dataset('aod_interpolated.nc', 'w', format='NETCDF4') as nc_out:
    nc_out.createDimension('points', len(longitude))
    lon_out = nc_out.createVariable('longitude', np.float32, ('points',))
    lon_out[:] = longitude
    lat_out = nc_out.createVariable('latitude', np.float32, ('points',))
    lat_out[:] = latitude
    aod_out = nc_out.createVariable('interpolated_aod', np.float32, ('points',))
    aod_out[:] = interpolated_aod
    obs_out = nc_out.createVariable('aeronet_aod', np.float32, ('points',))
    obs_out[:] = obs_aod

# Create a Lambert Conformal projection for plotting
lambert_proj = ccrs.LambertConformal(
    central_latitude=cen_lat,
    central_longitude=cen_lon,
    standard_parallels=(lat1, lat2)
)

# Plot the interpolated AOD
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

# Initialize arrays for transformed coordinates
longitude_lcc = np.zeros(longitude.shape)
latitude_lcc = np.zeros(latitude.shape)

# Transform the data point to projection coordinates
# Transform each point
for i in range(len(longitude)):
    longitude_lcc[i], latitude_lcc[i] = lambert_proj.transform_point(longitude[i], latitude[i], ccrs.PlateCarree())

scatter = ax.scatter(longitude_lcc, latitude_lcc, c=interpolated_aod, cmap='YlOrRd', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='Interpolated AOD')
plt.title('AOD Interpolated Values')
plt.savefig('aod_cmaq_interpolated.png')
plt.show()

# Plot the AERONET AOD
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

scatter = ax.scatter(longitude_lcc, latitude_lcc, c=obs_aod, cmap='YlOrRd', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='AERONET AOD')
plt.title('AOD AERONET Values')
plt.savefig('aod_aeronet_withnon.png')
plt.show()

# Plot the difference
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

scatter = ax.scatter(longitude_lcc, latitude_lcc, c=aod_difference, cmap='YlOrRd', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='Difference (Interpolated AOD - AERONET AOD)')
plt.title('Difference between Interpolated AOD and AERONET AOD')
plt.savefig('aod_difference.png')
plt.show()


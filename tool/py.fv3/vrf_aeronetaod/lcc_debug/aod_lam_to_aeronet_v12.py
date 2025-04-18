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
dx = 13000
dy = 13000
nx = 393
ny = 225

# Set up the Lambert projection with ellipsoid
lambert_proj = Proj(proj='lcc', lat_1=lat1, lat_2=lat2, lat_0=cen_lat, lon_0=cen_lon,
                     x_0=0, y_0=0, no_defs=True, ellps='WGS84')

# (0,0) is center in the above lambert_proj, where (0,0) is the left bottom corner of 2d aod
# need adjust the center of the 2d grid to (0,0) 
cen_x_adj = 0.5*(nx-1) 
cen_y_adj = 0.5*(ny-1) 

# Load CMAQAOD data in phy netcdf
aod_file = 'file1.nc'
with nc.Dataset(aod_file, 'r') as aod_nc:
    grid_xt = aod_nc.variables['grid_xt'][:]
    grid_yt = aod_nc.variables['grid_yt'][:]

# Distance (Meter) to grid number
    grid_xt =  grid_xt/dx
    grid_yt =  grid_yt/dy
    lon = np.zeros((ny, nx))
    lat = np.zeros((ny, nx))
    xt2d = np.zeros((ny, nx))
    yt2d = np.zeros((ny, nx))
    for i in range(ny):
        for j in range(nx):
            lon[i, j], lat[i, j] = lambert_proj(grid_xt[j], grid_yt[i], inverse=True)
#xt2d and yt2d are for interpolation use coordinate grid on LCC
            xt2d[i, j] = grid_xt[j]
            yt2d[i, j] = grid_yt[i] 
    aod_values = aod_nc.variables['aod'][0, :, :]

# Load AERONET data
aeronet_ds = xr.open_dataset('file2.nc')
obs_aod = aeronet_ds['aod_550nm'][0, 0].values
latitude = aeronet_ds['latitude'].values
longitude = aeronet_ds['longitude'].values
obs_aod[obs_aod == -1] = np.nan  # Replace missing values with NaN

#longitude = np.where(longitude < 0, longitude + 360, longitude)


#Extreme case test
#longitude[:] = -100.0
#latitude[:] = 38.5

print("min and max (longitude):")
#print(longitude)
print(longitude.min())
print(longitude.max())

# Print shape and size of obs_aod
print(f"obs_aod shape: {obs_aod.shape}, size: {obs_aod.size}")
print(f"latitude shape: {latitude.shape}, size: {latitude.size}")
#print("obs_aod:")
#print(aod_values)

#Test transforms between LCC grid and Lon,Lat coordinate 
#test_lon = cen_lon+1.0  # Example longitude
#test_lat = cen_lat+1.0  # Example latitude
#obs_x, obs_y = lambert_proj(test_lon, test_lat)
#obs_x = obs_x/dx
#obs_y = obs_y/dy
#print("Test Observation X coordinates (obs_x):", obs_x)
#print("Test Observation Y coordinates (obs_y):", obs_y)
#test_lon, test_lat =  lambert_proj(obs_x, obs_y, inverse=True)
#print("Test Observation X coordinates (obs_lon):", test_lon)
#print("Test Observation Y coordinates (obs_lat):", test_lat)

# Map observation coordinates to LCC coordinate grid
obs_x, obs_y = lambert_proj(longitude, latitude)
obs_x = obs_x/dx
obs_y = obs_y/dy

# Convert to cmaq model grid. Left bottom is (0,0)
obs_x = obs_x+cen_x_adj
obs_y = obs_y+cen_y_adj

print(f"obs_x shape: {obs_x.shape}, size: {obs_x.size}")
# Print obs_x and obs_y values
print("Observation X coordinates (obs_x):")
#print(obs_x)
print(obs_x.min())
print(obs_x.max())
#print("Observation Y coordinates (obs_y):")
#print(obs_y)
print(obs_y.min())
print(obs_y.max())

num_elements = len(obs_x)  # or len(obs_y), assuming both are the same length
valid_mask = np.full(num_elements, False)  # Create an array of True values

# Create a mask for valid observations
valid_mask = (obs_x >= grid_xt.min()+0.5) & (obs_x <= grid_xt.max()-0.5) & \
             (obs_y >= grid_yt.min()+0.5) & (obs_y <= grid_yt.max()-0.5)    

print("min,max (grid_xt,yt):")
print(grid_xt.min())
print(grid_xt.max())
print(grid_yt.min())
print(grid_yt.max())
# Print shape and size of valid_mask
print(f"valid_mask shape: {valid_mask.shape}, size: {valid_mask.size}")

# Count the number of True values in valid_mask
num_valid_obs = np.sum(valid_mask)

# Print the number of valid observations
print(f"Number of valid observations: {num_valid_obs}")

# Filter valid observations
valid_longitude = longitude[valid_mask]
valid_latitude = latitude[valid_mask]
valid_obs_aod = obs_aod[valid_mask]

# Print the valid longitude and latitude values
#print("Valid Longitude:")
#print(valid_longitude)

#print("Valid Latitude:")
#print(valid_latitude)

print("min and max (valid_longitude):")
print(valid_longitude.min())
print(valid_longitude.max())
print(valid_latitude.min())
print(valid_latitude.max())

# Create an array filled with NaNs
cmaq_aod_at_obsloc = np.full_like(obs_aod, np.nan)

# Interpolation using valid observations
for i in range(len(longitude)):
    print(f"lon= {longitude[i]} lat= {latitude[i]} ")
    print(f"obs= {obs_aod[i]} obs_x= {obs_x[i]} obs_y= {obs_y[i]} ")
    print("valid_mask= ", valid_mask[i])

    if valid_mask[i]:  # Use valid_mask[i] instead of valid_mask(i)
        ii = int(obs_x[i])
        jj = int(obs_y[i])
        
        print(f"X/Y location is: {obs_x[i]} {obs_y[i]}")
        print(f"Integer location is: {ii} {jj}")
        
        dx = obs_x[i] - ii
        dy = obs_y[i] - jj
        dx1 = 1.0 - dx
        dy1 = 1.0 - dy
        
        print(f"Offsets: {dx1} {dy1}")
        
        # Bilinear interpolation calculation
        cmaq_val = (
            aod_values[jj, ii] * dx1 * dy1 +
            aod_values[jj + 1, ii + 1] * dx * dy +
            aod_values[jj, ii + 1] * dx * dy1 +
            aod_values[jj + 1, ii] * dx1 * dy
        )
        
        print(f"obs= {obs_aod[i]} cmaq_val= {cmaq_val} "
              f"longitude= {longitude[i]} latitude= {latitude[i]} "
              f"ii= {ii} jj= {jj} "
              f"aod_values[jj, ii]= {aod_values[jj, ii]} "
              f"aod_values[jj + 1, ii + 1]= {aod_values[jj + 1, ii + 1]}")
        
        cmaq_aod_at_obsloc[i] = cmaq_val

# Interpolation using valid observations
interpolated_aod_a = np.full_like(obs_aod, np.nan)
if valid_longitude.size > 0:
    interpolated_aod_a[valid_mask] = griddata(
        (lon.flatten(), lat.flatten()),
        aod_values.flatten(),
        (valid_longitude, valid_latitude),
        method='linear',
        fill_value=np.nan
    )

interpolated_aod = np.full_like(obs_aod, np.nan)
if valid_longitude.size > 0:
    interpolated_aod = griddata(
        (xt2d.flatten(), yt2d.flatten()),
        aod_values.flatten(),
        (obs_x, obs_y),
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
print(f"Mean: {np.nanmean(valid_obs_aod):.3f}, Std: {np.nanstd(valid_obs_aod):.3f}, "
      f"Min: {np.nanmin(valid_obs_aod):.3f}, Max: {np.nanmax(valid_obs_aod):.3f}")

print("Interpolated AOD_LCC Statistics:")
print(f"Mean: {np.nanmean(interpolated_aod):.3f}, Std: {np.nanstd(interpolated_aod):.3f}, "
      f"Min: {np.nanmin(interpolated_aod):.3f}, Max: {np.nanmax(interpolated_aod):.3f}")

print("Interpolated AOD_LL Statistics:")
print(f"Mean: {np.nanmean(interpolated_aod_a):.3f}, Std: {np.nanstd(interpolated_aod_a):.3f}, "
      f"Min: {np.nanmin(interpolated_aod_a):.3f}, Max: {np.nanmax(interpolated_aod_a):.3f}")

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

scatter = ax.scatter(longitude_lcc, latitude_lcc, c=interpolated_aod, cmap='coolwarm', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='Interpolated AOD')
plt.title('AOD Interpolated Values')
plt.savefig('aod_interpolated.png')
plt.show()

# Plot the AERONET AOD
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

#scatter = ax.scatter(valid_longitude, valid_latitude, c=valid_obs_aod, cmap='coolwarm', edgecolor='k', alpha=0.7)
scatter = ax.scatter(longitude_lcc, latitude_lcc, c=obs_aod, cmap='coolwarm', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='AERONET AOD')
plt.title('AOD AERONET Values')
plt.savefig('aod_aeronet.png')
plt.show()

# Plot the difference
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

scatter = ax.scatter(longitude_lcc, latitude_lcc, c=aod_difference, cmap='coolwarm', edgecolor='k', alpha=0.7)
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')
#ax.gridlines(draw_labels=True)
plt.colorbar(scatter, label='Difference (Interpolated AOD - AERONET AOD)')
plt.title('Difference between Interpolated AOD and AERONET AOD')
plt.savefig('aod_difference.png')
plt.show()


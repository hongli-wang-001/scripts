import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.crs as ccrs

# Open NetCDF file
nc_path = 'combined_ioda.nc'
ds = Dataset(nc_path, mode='r')

# Access the MetaData group
meta = ds.groups['MetaData']
obs = ds.groups['ObsValue']

# Read variables
lat = meta.variables['latitude'][:]
lon = meta.variables['longitude'][:]
zenith = meta.variables['sensorZenithAngle'][:]
azimuth = meta.variables['sensorAzimuthAngle'][:]
cloud = meta.variables['cloudAmount'][:]

rec_rad = obs.variables['radiance'][:,:]
# Mask invalid values
valid_mask = (rec_rad >= 0.0) 
# Use masked array so that invalid values are ignored
masked = np.ma.array(rec_rad, mask=~valid_mask)

# Compute stats along the third dimension (axis = 2)
mean_per_channel = masked.mean(axis=(0))       # shape → (dim3,)
median_per_channel = np.ma.median(masked, axis=(0))
std_per_channel = masked.std(axis=(0))         # default ddof=0; use ddof=1 for sample std

# Print / inspect
for i in range(rec_rad.shape[1]):
    print(f"Channel {i}: mean={mean_per_channel[i]:.4f}, "
          f"median={median_per_channel[i]:.4f}, std={std_per_channel[i]:.4f}")

radiance = obs.variables['radiance'][:,149]
range_mask = (radiance < 0.0) | (radiance > 20.0)
radiance = np.ma.masked_where(range_mask, radiance)

# Handle fill values
def mask_fill(var):
    fill = var._FillValue if hasattr(var, '_FillValue') else None
    data = var[:]
    if fill is not None:
        data = np.ma.masked_where(data == fill, data)
    return data

# Read and mask lat/lon
lat = mask_fill(meta.variables['latitude'])
lon = mask_fill(meta.variables['longitude'])

# Mask invalid placeholder values (0.6554E+05)
invalid_val = 65540.0
mask_invalid = (lat == invalid_val) | (lon == invalid_val)
lat = np.ma.masked_where(mask_invalid, lat)
lon = np.ma.masked_where(mask_invalid, lon)

# Read and mask sensor zenith angle
zenith = mask_fill(meta.variables['sensorZenithAngle'])
zenith = np.ma.masked_where(mask_invalid, zenith)  # apply same mask

range_mask = (radiance < 0.0) | (radiance > 20.0)
total_mask = mask_invalid | range_mask

radiance = np.ma.masked_where(total_mask, radiance)

#radiance = np.ma.masked_where(mask_invalid, radiance)  # apply same mask

azimuth = mask_fill(meta.variables['sensorAzimuthAngle'])
cloud = mask_fill(meta.variables['cloudAmount'])


#Lat/lon range
latlonrange = np.fromstring("0 0 0 0", dtype=float, sep=' ')
latlonrange[0] = np.amin(lon)
latlonrange[1] = np.amax(lon)
latlonrange[2] = np.amin(lat)
latlonrange[3] = np.amax(lat)

# Quick scatter plot (color-coded by zenith angle)
plt.figure(figsize=(12, 9))
ax = plt.axes(projection=ccrs.PlateCarree())

# Set extent and map features
ax.set_extent([-90, 90, 0, 90], crs=ccrs.PlateCarree())
ax.coastlines('50m')

sc = plt.scatter(lon, lat, c=zenith, cmap='viridis', s=1, alpha=0.5)
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='gray', linestyle='--', draw_labels=True)
gl.top_labels = None
# Plot scatter on the Cartopy axis
# cmap='viridis'
cmap="jet"       #alternative: "RdBu_r"
sc = ax.scatter(lon, lat, c=zenith, cmap=cmap, s=1, alpha=0.5, transform=ccrs.PlateCarree())
plt.colorbar(sc, label='Sensor Zenith Angle (°)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Observation Locations Colored by Sensor Zenith Angle')
plt.xlim(-90,90)
plt.ylim(0,90)
plt.grid(True)
plt.tight_layout()
#plt.show()
# Save the figure
outfig="sensorzenith.png"
plt.savefig(outfig, dpi=300)
exit

# Create figure and Cartopy axis
fig = plt.figure(figsize=(12, 9))
ax = plt.axes(projection=ccrs.PlateCarree())

# Set extent and map features
ax.set_extent([-90, 90, 0, 90], crs=ccrs.PlateCarree())
ax.coastlines('50m')
gl = ax.gridlines(crs=ccrs.PlateCarree(), linewidth=1, color='gray', linestyle='--', draw_labels=True)
gl.top_labels = None

# Plot scatter on the Cartopy axis
# cmap='viridis'
cmap="jet"       #alternative: "RdBu_r"
sc = ax.scatter(lon, lat, c=radiance, cmap=cmap, s=1, alpha=0.5, transform=ccrs.PlateCarree())

#plt.scatter(lon,lat,s=int(args.symsize),c=rad, marker="o", lw=0, cmap=cmap)

# Draw the colour bar
#cbar = plt.colorbar(orientation='horizontal', shrink = 0.6)
#cbar.ax.set_xlabel('radiance (W/m$^2$/sr/cm$^{-1}$)')



plt.colorbar(sc, label='Radiance')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Reconstructed Radiance (unit: W m-2 sr-1)')
plt.xlim(-90,90)
plt.ylim(0,90)
plt.grid(True)
plt.tight_layout()
#plt.show()
# Save the figure
outfig="radiance_v3.png"
plt.savefig(outfig, dpi=300)

# Close the dataset
ds.close()


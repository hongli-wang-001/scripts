import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from netCDF4 import Dataset
import re
from datetime import datetime

states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
# -------------------------------
# Set file path and extract time
# -------------------------------
filename = 'pm25.2020091000_f024.nc'  # Replace with your actual file
match = re.search(r'\.(\d{10})', filename)
if match:
    time_str = match.group(1)  # e.g., '2020091000'
    dt = datetime.strptime(time_str, '%Y%m%d%H')
    title_time = dt.strftime('%Y-%m-%d %H:00 UTC')
else:
    title_time = "Unknown time"

# -------------------------------
# Load data
# -------------------------------
ds = Dataset(filename)
lat = ds.variables['latitude'][:]     # shape (225, 393)
lon = ds.variables['longitude'][:]    # shape (225, 393)
pm25 = ds.variables['pm25'][:]        # shape (225, 393)
pm25_aqm = ds.variables['pm25_aqm'][:]   

# -------------------------------
# Define contour levels and colormap
# -------------------------------
levels = [0.1, 2, 3, 5, 7, 10, 15, 20, 35, 50, 70, 100]
cmap = 'jet'

# -------------------------------
# Plotting
# -------------------------------
plt.figure(figsize=(12, 6))
ax=plt.axes(projection=ccrs.PlateCarree())
#cf = plt.contourf(lon, lat, pm25, levels=levels, cmap=cmap, extend='max')
cf = plt.scatter(
    lon.flatten(), lat.flatten(), s=2, c=pm25.flatten(), cmap='jet', 
    norm=colors.BoundaryNorm(boundaries=levels, ncolors=256)
)
#plt.axis([-135,-58,24,55])
cb = plt.colorbar(cf, label='PM2.5 (µg/m³)')
cb.set_ticks(levels)
cb.set_ticklabels([str(l) for l in levels])
ax.set_extent([225, 300, 24, 55], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title(f'AirNow PM2.5\n{title_time}')
plt.savefig('obs_pm25_h2d.png',dpi=300)
plt.tight_layout()
plt.show()


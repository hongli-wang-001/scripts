import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from netCDF4 import Dataset
import re
from datetime import datetime

# -------------------------------
# Define state/province borders
# -------------------------------
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none'
)

# -------------------------------
# Set file path and extract time
# -------------------------------
filename = 'pm25.2020090201_f001.nc'  # Replace with your actual file
match = re.search(r'\.(\d{10})', filename)
if match:
    time_str = match.group(1)
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
dsb = Dataset('i18/statistics.f001.nc')
pm25 = dsb.variables['mean_pm25_con'][0,:,:]        # observations
pm25_aqm = dsb.variables['mean_pm25_con'][0,:,:]  # model

# Flatten arrays
flat_lon = lon.flatten()
flat_lat = lat.flatten()
flat_pm25 = pm25.flatten()
# Mask out NaNs
valid_mask = ~np.isnan(flat_pm25)

# Extract Lambert Conformal projection parameters
proj_params = ds.variables['Lambert_Conformal'] if 'Lambert_Conformal' in ds.variables else ds

#    cen_lat	: 	38.5
#    cen_lon	: 	-97.5
#    dx	: 	13000
#    dy	: 	13000
#    grid	: 	lambert_conformal
#    stdlat1	: 	38.5
#    stdlat2	: 	38.5

std_parallels = proj_params.standard_parallel[:] if hasattr(proj_params, 'standard_parallel') else [38.5, 38.5]
central_lon = proj_params.longitude_of_central_meridian if hasattr(proj_params, 'longitude_of_central_meridian') else -97.5
lat_origin = proj_params.latitude_of_projection_origin if hasattr(proj_params, 'latitude_of_projection_origin') else 38.5

# -------------------------------
# Define AQI-style levels and colormap
# -------------------------------
#aqi_levels = [0.0, 12.1, 35.5, 55.5, 150.5, 250.5, 500.0]
#aqi_colors = ['#00e400', '#ffff00', '#ff7e00', '#ff0000', '#8f3f97', '#7e0023']
aqi_colors = ['#00e400', '#ffff00', '#ff7e00', '#ff0000', '#8f3f97', '#7e0023']
# New levels with intermediate ones
#aqi_levels = [0.0, 0.1, 1.0, 2.0, 4, 6, 8, 12.1, 35.5, 55.5, 150.5, 250.5, 500.0]
aqi_levels = [100, 200.0, 300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1250., 1500.0, 2000.0, 2500.0]
# New color mapping: the lower values have different shades, but above 12.1 it stays the same.
aqi_colors = [
#    '#00e400',  # 0.0–0.1   (Green - Good)
#    '#7ee600',  # 0.1–1.0   (Light green)
#    '#a8e05f',  # 1.0–2.0   (Lime)
#    '#d6ef39',  # 2.0–4.0   (Yellow-green)
#    '#ffff00',  # 4.0–6.0   (Yellow - Moderate)
#    '#ffcc00',  # 6.0–8.0   (Yellow-orange)
#    '#ff9933',  # 8.0–12.1  (Orange - Unhealthy for Sensitive Groups)
    'white',    # 0.0–1.0
    '#d0f0ff',  # 1.0–2.0   very light blue
    '#a0e8ff',  # 2.1–4.0   light blue
    '#70d8ff',  # 4.0–6.0   sky blue
    '#40c0ff',  # 6.0–8.0   cyan-blue

    '#00e400',  # 8.0–12.1  green (standard "Good" color
    '#ffff00',  # 12.1–35.5 (Yellow - Unhealthy)
    '#ff7e00',  # 35.5–55.5 (Orange - Unhealthy)
    '#ff0000',  # 55.5–150.5 (Red - Very Unhealthy)
    '#8f3f97',  # 150.5–250.5 (Purple - Very Unhealthy)
    '#7e0023',  # 250.5–500.0 (Maroon - Hazardous)
]

aqi_cmap = mcolors.ListedColormap(aqi_colors)
aqi_norm = mcolors.BoundaryNorm(aqi_levels, len(aqi_colors))

# -------------------------------
# Plotting
# -------------------------------
proj = ccrs.LambertConformal(central_longitude=central_lon, standard_parallels=std_parallels, 
                             central_latitude=lat_origin)
fig = plt.figure(figsize=(12, 6))
ax = plt.axes(projection=proj)

# Filled contour for model
contour = ax.contourf(
    lon, lat, pm25_aqm, levels=aqi_levels, cmap=aqi_cmap, norm=aqi_norm,
    transform=ccrs.PlateCarree(), extend='max'
)

# Overlay observations as scatter
#scatter = ax.scatter(
#    flat_lon[valid_mask], flat_lat[valid_mask], s=10, c=flat_pm25[valid_mask],
#    cmap=aqi_cmap, norm=aqi_norm, edgecolor='k', linewidth=0.3,
##    edgecolors='none',  # avoids black borders around points
#    transform=ccrs.PlateCarree()
#)

# Map features
#ax.set_extent([-125, -66.5, 24, 50], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')

# Colorbar
cb = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.02, shrink=0.85)
cb.set_label('PBLH (m)')
cb.set_ticks(aqi_levels)
#cb.set_ticklabels(['Good', 'Moderate', 'USG', 'Unhealthy', 'Very Unhealthy', 'Hazardous', 'Hazardous'])

# Title and layout
plt.title(f'PHLB')
plt.tight_layout()
plt.savefig('pblh_i18.png', dpi=300)
plt.show()


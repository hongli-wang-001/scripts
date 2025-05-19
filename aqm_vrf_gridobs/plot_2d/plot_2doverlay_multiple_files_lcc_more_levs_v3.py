import os
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from netCDF4 import Dataset
import re
from datetime import datetime

states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none'
)

aqi_levels = [0.0, 1.0, 2.0, 4, 6, 8, 12.1, 35.5, 55.5, 150.5, 250.5, 500.0]
aqi_colors = [
    'white', '#d0f0ff', '#a0e8ff', '#70d8ff', '#40c0ff',
    '#00e400', '#ffff00', '#ff7e00', '#ff0000', '#8f3f97', '#7e0023'
]
aqi_cmap = mcolors.ListedColormap(aqi_colors)
aqi_norm = mcolors.BoundaryNorm(aqi_levels, len(aqi_colors))

file_list = sorted(glob.glob("pm25.20*.nc"))
for filename in file_list:
    match = re.search(r'\.(\d{10})', filename)
    if match:
        time_str = match.group(1)
        dt = datetime.strptime(time_str, '%Y%m%d%H')
        title_time = dt.strftime('%Y-%m-%d %H:00 UTC')
    else:
        title_time = "Unknown time"

    ds = Dataset(filename)
    lat = ds.variables['latitude'][:]
    lon = ds.variables['longitude'][:]
    pm25 = ds.variables['pm25'][:]
    pm25_aqm = ds.variables['pm25_aqm'][:]

    flat_lon = lon.flatten()
    flat_lat = lat.flatten()
    flat_pm25 = pm25.flatten()
    valid_mask = ~np.isnan(flat_pm25)

    proj_params = ds.variables['Lambert_Conformal'] if 'Lambert_Conformal' in ds.variables else ds
    std_parallels = proj_params.standard_parallel[:] if hasattr(proj_params, 'standard_parallel') else [38.5, 38.5]
    central_lon = proj_params.longitude_of_central_meridian if hasattr(proj_params, 'longitude_of_central_meridian') else -97.5
    lat_origin = proj_params.latitude_of_projection_origin if hasattr(proj_params, 'latitude_of_projection_origin') else 38.5

    proj = ccrs.LambertConformal(central_longitude=central_lon, standard_parallels=std_parallels,
                                 central_latitude=lat_origin)

    fig = plt.figure(figsize=(12, 6))
    ax = plt.axes(projection=proj)

    contour = ax.contourf(
        lon, lat, pm25_aqm, levels=aqi_levels, cmap=aqi_cmap, norm=aqi_norm,
        transform=ccrs.PlateCarree(), extend='max'
    )
    scatter = ax.scatter(
        flat_lon[valid_mask], flat_lat[valid_mask], s=10, c=flat_pm25[valid_mask],
        cmap=aqi_cmap, norm=aqi_norm, edgecolor='k', linewidth=0.3,
        transform=ccrs.PlateCarree()
    )

    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.COASTLINE)
    ax.add_feature(states_provinces, edgecolor='gray')

    cb = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.02, shrink=0.85)
    cb.set_label('PM2.5 (µg/m³)')
    cb.set_ticks(aqi_levels)

    plt.title(f'PM2.5: Model (Shaded) + Observations @ {title_time}')
    plt.tight_layout()
    plt.savefig(f"obs_vs_model_pm25_{time_str}.png", dpi=300)
    plt.close()


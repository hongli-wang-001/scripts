import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature

# Define Lambert conformal projection parameters
lat1 = 22.57418
lat2 = 38.5
cen_lat = 38.5
cen_lon = -97.5

# Create a Lambert Conformal projection for plotting
lambert_proj = ccrs.LambertConformal(
    central_latitude=cen_lat,
    central_longitude=cen_lon,
    standard_parallels=(lat1, lat2)
)

# Create a plot
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())

# Add features to the map
ax.coastlines()
ax.add_feature(cartopy.feature.BORDERS)
ax.add_feature(cartopy.feature.STATES)

# Example data point
lon_point = -100  # Example longitude
lat_point = 35    # Example latitude

# Transform the data point to projection coordinates
obs_x, obs_y = lambert_proj.transform_point(lon_point, lat_point, ccrs.PlateCarree())

# Plot the data point
ax.scatter(obs_x, obs_y, color='red', s=100, marker='o', label='Data Point')

plt.title("Lambert Conformal Projection")
plt.legend()
plt.savefig('oneobs_on_lcc.png')

plt.show()


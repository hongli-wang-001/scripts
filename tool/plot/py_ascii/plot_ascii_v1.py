import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Load the data from the ASCII file
data = np.loadtxt('data.txt')  # Replace 'data.txt' with your actual file name

# Extract longitude, latitude, and AOD (4th column)
longitudes = data[:, 0]  # 1st column
latitudes = data[:, 1]   # 2nd column
aod_values = data[:, 3]  # 4th column

# Create a Lambert Conformal projection for plotting
cen_lat = 38.5
cen_lon = -97.5
lambert_proj = ccrs.LambertConformal(central_latitude=cen_lat, central_longitude=cen_lon)

# Create a plot
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())  # Set the extent of the map

# Add map features
ax.add_feature(cfeature.COASTLINE, linestyle='--', color='blue')  # Add coastlines
ax.add_feature(cfeature.BORDERS, linestyle=':', color='black')    # Add borders
ax.add_feature(cfeature.STATES, linestyle='-', color='black')     # Add state borders

# Add gridlines
ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')

# Transform the coordinates to Lambert projection
obs_x, obs_y = lambert_proj.transform_points(ccrs.PlateCarree(), longitudes, latitudes).T

# Plot the AOD data
scatter = ax.scatter(obs_x, obs_y, c=aod_values, cmap='coolwarm', edgecolor='k', alpha=0.7, s=50)

# Add a colorbar
cbar = plt.colorbar(scatter, ax=ax, orientation='vertical', pad=0.01)
cbar.set_label('AOD Values')

plt.title("AOD Values on Lambert Conformal Projection")
plt.savefig('aod_map.png')  # Save the figure
plt.show()


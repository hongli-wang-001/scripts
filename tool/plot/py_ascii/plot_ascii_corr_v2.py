import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Load the data from the ASCII file
data = np.loadtxt('data.txt')  # Replace 'data.txt' with your actual file name

# Extract longitude, latitude, and bias (4th column, assuming it's bias here)
longitudes = data[:, 0]  # 1st column
latitudes = data[:, 1]   # 2nd column
bias_values = data[:, 6]  # 4th column (assuming this is the bias)
bias_values[bias_values == -99.0] = np.nan 

# Create a Lambert Conformal projection for plotting
cen_lat = 38.5
cen_lon = -97.5
lambert_proj = ccrs.LambertConformal(central_latitude=cen_lat, central_longitude=cen_lon)

# Create a plot
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=lambert_proj)
ax.set_extent([-125, -70, 25, 50], ccrs.PlateCarree())  # Set the extent of the map

# Add map features
ax.add_feature(cfeature.COASTLINE, linestyle='--')
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.STATES, linestyle='-')

# Add gridlines
ax.gridlines(draw_labels=True, color='gray', alpha=0.5, linestyle='--')

# Transform the coordinates to Lambert projection
transformed_points = lambert_proj.transform_points(ccrs.PlateCarree(), longitudes, latitudes)
obs_x = transformed_points[:, 0]
obs_y = transformed_points[:, 1]

# Create a symmetric color map centered around white
vmin = bias_values.min()
vmax = bias_values.max()
norm = plt.Normalize(vmin=-max(abs(vmin), abs(vmax)), vmax=max(abs(vmin), abs(vmax)))

# Plot the bias data
scatter = ax.scatter(obs_x, obs_y, c=bias_values, cmap='RdBu_r', norm=norm, edgecolor='k', alpha=0.7, s=50)

# Add a colorbar
cbar = plt.colorbar(scatter, ax=ax, orientation='vertical', pad=0.06,shrink=0.7)
cbar.set_label('Bias Values')

plt.title("Correlation Values")
plt.savefig('corr_map.png')  # Save the figure
plt.show()


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from netCDF4 import Dataset
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def plot_particulate_matter(nc_file, output_file='pm25_surface_plot.png'):
    # Open the NetCDF file
    with Dataset(nc_file, 'r') as nc:
        # Access latitude and longitude from the MetaData group
        lat = nc.groups['MetaData'].variables['latitude'][:]
        lon = nc.groups['MetaData'].variables['longitude'][:]

        # Access particulate matter from the ObsValue group
        pm25_surface = nc.groups['ObsValue'].variables['particulatematter2p5Surface'][:]

        # Handle fill values
        fill_value = nc.groups['ObsValue'].variables['particulatematter2p5Surface']._FillValue
        pm25_surface = np.where(pm25_surface == fill_value, np.nan, pm25_surface)

    # Create a map with Cartopy
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([-140, -58, 24, 55], ccrs.PlateCarree())
    
    # Add geographical features
    ax.add_feature(cfeature.COASTLINE, linestyle='--')
    ax.add_feature(cfeature.BORDERS, linestyle=':')
    ax.add_feature(cfeature.STATES, linestyle='-')

    # Scatter plot
    myvalues = [2, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100]
    #scatter = ax.scatter(lon, lat, c=pm25_surface, cmap='viridis', marker='o', edgecolor='k', s=20)
    #sc = plt.scatter(tlon, tlat, s=5, c=pm, linewidth=.50, cmap='jet', norm=colors.BoundaryNorm(boundaries=myvalues, ncolors=256))
    scatter = ax.scatter(lon, lat, s=3, c=pm25_surface, linewidth=.50, cmap='jet', norm=colors.BoundaryNorm(boundaries=myvalues, ncolors=256))
    # Colorbar and labels
    plt.colorbar(scatter, label='PM2.5 Surface (µg m⁻³)')
    plt.title('PurpleAir PM2.5')
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('Latitude (degrees)')
    plt.grid()

    # Save the figure
    plt.savefig(output_file)
    plt.show()

# Replace 'your_file.nc' with the path to your NetCDF file
plot_particulate_matter('file1.nc', output_file='pm25_surface_plot.png')


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

def plot_particulate_matter(nc_file, output_file='pm25_surface_plot.png'):
    # Open the NetCDF file
    with Dataset(nc_file, 'r') as nc:
        # Read longitude, latitude, and particulatematter2p5Surface
        longitude = nc.variables['longitude'][:]
        latitude = nc.variables['latitude'][:]
        pm25_surface = nc.variables['particulatematter2p5Surface'][:]

        # Handle fill values
        fill_value = nc.variables['particulatematter2p5Surface']._FillValue
        pm25_surface = np.where(pm25_surface == fill_value, np.nan, pm25_surface)

    # Create a scatter plot
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(longitude, latitude, c=pm25_surface, cmap='viridis', marker='o', edgecolor='k', s=20)
    plt.colorbar(scatter, label='PM2.5 Surface (µg m⁻³)')
    plt.title('Particulate Matter (PM2.5) Surface Concentration')
    plt.xlabel('Longitude (degrees)')
    plt.ylabel('Latitude (degrees)')
    plt.grid()

    # Save the figure
    plt.savefig(output_file)
    plt.show()

# Replace 'your_file.nc' with the path to your NetCDF file
plot_particulate_matter('file1.nc', output_file='pm25_surface_plot.png')


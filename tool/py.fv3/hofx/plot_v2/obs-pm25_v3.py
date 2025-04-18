#!/usr/bin/env python
# Plot PM2.5 with custom colormap

import netCDF4 as nc
import sys
import numpy as np  # Import numpy for calculations
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Define the custom colormap
def create_custom_cmap():
    """Create a custom colormap."""
    cmap = mcolors.ListedColormap([
        (0.0000, 0.7060, 0.0000), (0.0000, 0.9060, 0.0000), (0.3020, 1.0000, 0.3020),
        (1.0000, 1.0000, 0.4980), (1.0000, 0.8745, 0.0000), (1.0000, 0.6471, 0.0000),
        (1.0000, 0.3840, 0.3840), (1.0000, 0.0000, 0.0000), (0.8000, 0.0000, 0.0000), 
        (0.7020, 0.0000, 0.0000), (0.6120, 0.5100, 0.8120), (0.5180, 0.3880, 0.7650),
        (0.4310, 0.2780, 0.7250), (0.2980, 0.1920, 0.5020), (0.4706, 0.4706, 0.4706),
        (0.7843, 0.7843, 0.7843)
    ])
    cmap.set_under((0.8627, 0.8627, 1.0000))  # Color for values below min
    cmap.set_over((0.9412, 0.9412, 0.9412))   # Color for values above max
    return cmap

# Ensure usage
if len(sys.argv) < 4:
    print('Wrong usage:')
    print(sys.argv[0] + ' infile title outfig')
    sys.exit(1)

# Read data
f1 = nc.Dataset(sys.argv[1], 'r')
tlon = f1['MetaData/longitude'][:]
tlat = f1['MetaData/latitude'][:]
pm = f1['ObsValue/particulatematter2p5Surface'][:]

# Calculate statistics
min_pm = np.min(pm)
max_pm = np.max(pm)
mean_pm = np.mean(pm)
median_pm = np.median(pm)
std_pm = np.std(pm)

myvalues = [ 3., 6., 9., 12., 15., 35., 55., 75., 100., 125., 150., 250., 300., 400., 500., 600., 750. ]  
# Create the plot
fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': ccrs.PlateCarree()})
cmap = create_custom_cmap()

# Scatter plot
scatter = ax.scatter(tlon, tlat, s=5, c=pm, linewidth=0.50, cmap=cmap, edgecolor='none', norm=mcolors.BoundaryNorm(boundaries=myvalues, ncolors=256))

# Add features
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none'), edgecolor='gray')

# Add colorbar
cbar = plt.colorbar(scatter, ax=ax, orientation='horizontal', fraction=0.05, pad=0.1)
cbar.set_label('Surface PM2.5 ($\mu$g/m$^3$)', fontsize=16)

# Add text with data statistics
stats_text = (f'Min: {min_pm:.2f}\nMax: {max_pm:.2f}\nMean: {mean_pm:.2f}\n'
              f'Median: {median_pm:.2f}\nStd Dev: {std_pm:.2f}')
plt.text(0.05, 0.05, stats_text, transform=ax.transAxes, fontsize=12,
         verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.5'))

# Plot title and axis labels
plt.title(sys.argv[2], fontsize=18)

# Set axis limits
plt.axis([-140, -58, 24, 55])

# Adjust layout and save figure
plt.tight_layout()
plt.savefig(sys.argv[3], dpi=300)
plt.show()


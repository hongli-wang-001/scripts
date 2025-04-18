#!/usr/bin/env python
# plot HofX 
import netCDF4 as nc
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Function to compute the Modified Z-Score
def modified_z_score(data):
    median = np.median(data)
    mad = np.median(np.abs(data - median))
    # Handle the case where MAD is zero (no variation)
    mad = mad if mad != 0 else 1
    return 0.6745 * (data - median) / mad

if len(sys.argv) < 4:
    print('Wrong usage:')
    print(sys.argv[0] + ' infile title outfig')
    sys.exit(1)

myvalues = [2, 3, 5, 7, 10, 15, 20, 30, 50, 70, 100]

# Read data from NetCDF file
f1 = nc.Dataset(sys.argv[1], 'r')
tlon = f1['MetaData/longitude'][:]
tlat = f1['MetaData/latitude'][:]
pm = f1['ObsValue/particulatematter2p5Surface'][:]

# Flatten the data for outlier detection
pm_flat = pm.flatten()

# Compute the Modified Z-Score
mzs = modified_z_score(pm_flat)

# Define a threshold for outlier detection
threshold = 3.5  # Commonly used threshold for outliers
outlier_mask = np.abs(mzs) > threshold

# Reshape the mask to the original shape
outlier_mask_reshaped = outlier_mask.reshape(pm.shape)

# Create a new figure for plotting
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

# Plot PM2.5 data
sc = plt.scatter(tlon, tlat, s=5, c=pm, linewidth=.50, cmap='jet', norm=colors.BoundaryNorm(boundaries=myvalues, ncolors=256))

# Highlight outliers
plt.scatter(tlon[outlier_mask], tlat[outlier_mask], s=20, c='red', marker='x', label='Outliers')

plt.axis([-140, -58, 24, 55])
cbar = plt.colorbar(sc, fraction=0.07, orientation='horizontal')
cbar.set_label('Surface PM2.5 ($\mu$g/m$^3$)', fontsize=16)
plt.title(sys.argv[2], fontsize=18)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none'
), edgecolor='gray')
plt.legend()
plt.tight_layout()
plt.savefig(sys.argv[3], dpi=300)


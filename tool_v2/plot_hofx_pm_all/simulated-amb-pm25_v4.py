import netCDF4 as nc
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

# Define states and provinces feature for Cartopy
states_provinces = cfeature.NaturalEarthFeature(
    category='cultural',
    name='admin_1_states_provinces_lines',
    scale='50m',
    facecolor='none'
)

# Define the levels for color mapping
myvalues = [-20, -10, -5, -3, -1, 1, 3, 5, 10, 20]

if len(sys.argv) < 4:
    print('Wrong usage:')
    print(sys.argv[0] + ' infile title outfig')
    sys.exit(1)

# Read the NetCDF file
infile = sys.argv[1]
title = sys.argv[2]
outfig = sys.argv[3]

try:
    with nc.Dataset(infile, 'r') as f1:
        # Extract longitude, latitude, and PM2.5 data
        tlon = f1['MetaData/longitude'][:]
        tlat = f1['MetaData/latitude'][:]
        pm = f1['hofx1/particulatematter2p5Surface'][:] - f1['hofx0/particulatematter2p5Surface'][:]
except Exception as e:
    print(f'Error reading the NetCDF file: {e}')
    sys.exit(1)

# Create the plot
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
scatter = plt.scatter(tlon, tlat, s=5, c=pm, linewidth=.50,
                      cmap='jet', norm=colors.BoundaryNorm(boundaries=myvalues, ncolors=256))

# Set axis limits
plt.axis([-140, -58, 24, 55])

# Add colorbar
cbar = plt.colorbar(scatter, fraction=0.07, orientation='horizontal')
cbar.set_label('Surface PM2.5 ($\mu$g/m$^3$)', fontsize=16)

# Add title and features
plt.title(title, fontsize=18)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')

# Add gridlines with labels
#gl = ax.gridlines(draw_labels=True, linestyle='--', color='gray', alpha=0.5)
gl = ax.gridlines(draw_labels=True, linestyle='--', color='none', alpha=0.0)
gl.xlabels_top = False
gl.ylabels_right = False
gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

# Save the figure
plt.tight_layout()
plt.savefig(outfig, dpi=300)

print(f'Plot saved as {outfig}')


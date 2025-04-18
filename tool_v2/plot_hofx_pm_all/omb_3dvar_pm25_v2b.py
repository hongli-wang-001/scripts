#!/usr/bin/env python
# plot HofX 
import netCDF4 as nc
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
	
if len(sys.argv) < 4:
    print('Wrong usage:')
    print(sys.argv[0] + ' infile title outfig')
    sys.exit(1)

# Statistics calculations
def print_stats(data, label):
    print(f"{label} - Mean: {np.mean(data):.2f}, Median: {np.median(data):.2f}, Std Dev: {np.std(data):.2f}, RMS: {np.sqrt(np.mean(data**2)):.2f}")

#myvalues=[2,3, 5, 7, 10, 15, 20, 30, 50, 70, 100]
myvalues=[-20, -10, -5, -3, -1, 1, 3, 5, 10, 20]

f1=nc.Dataset(sys.argv[1],'r')
tlon=f1['MetaData/longitude'][:]
tlat=f1['MetaData/latitude'][:]
pm=f1['ObsValue/particulatematter2p5Surface'][:]-f1['hofx0/particulatematter2p5Surface'][:]
print_stats(pm, "OmB")
ax=plt.axes(projection=ccrs.PlateCarree())
plt.scatter(tlon,tlat,s=5,c=pm,\
   linewidth=.50,cmap='jet',norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
plt.axis([-135,-58,24,55])   
cbar=plt.colorbar(fraction=0.07,orientation='horizontal')
cbar.set_label('Surface PM2.5 ($\mu$g/m$^3$)',fontsize=16)
plt.title(sys.argv[2],fontsize=18)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')
plt.tight_layout()
plt.savefig(sys.argv[3],dpi=300)

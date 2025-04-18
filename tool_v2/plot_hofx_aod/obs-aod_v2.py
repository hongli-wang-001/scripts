#!/usr/bin/env python
# plot HofX 
import netCDF4 as nc
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as colors
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
    print(f"{label}  - Min: {np.min(data):.2f},  - Max: {np.max(data):.2f}, - Mean: {np.mean(data):.2f}, Median: {np.median(data):.2f}, Std Dev: {np.std(data):.2f}, RMS: {np.    sqrt(np.mean(data**2)):.2f},  - Num: {np.size(data):.2f}")

myvalues=[0.01,0.01,0.05,0.1, 0.3,0.4,0.5, 0.7, 1.0, 1.5]
#myvalues=[-20, -10, -5, -3, -1, 1, 3, 5, 10, 20]

f1=nc.Dataset(sys.argv[1],'r')
tlon=f1['MetaData/longitude'][:]
tlat=f1['MetaData/latitude'][:]
#pm=f1['ObsValue/particulatematter2p5Surface'][:]
oaod=f1['ObsValue/aerosolOpticalDepth'][:,:].squeeze()
baod=f1['hofx0/aerosolOpticalDepth'][:,:].squeeze()
aaod=f1['hofx1/aerosolOpticalDepth'][:,:].squeeze()
pm=oaod
print_stats(pm, "OBS")

ax=plt.axes(projection=ccrs.PlateCarree())
plt.scatter(tlon,tlat,s=2,c=pm,\
   linewidth=.50,cmap='jet',norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
plt.axis([-135,-58,24,55])   
cbar=plt.colorbar(fraction=0.07,orientation='horizontal')
cbar.set_label('AOD',fontsize=16)
plt.title(sys.argv[2],fontsize=18)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')
plt.tight_layout()
plt.savefig(sys.argv[3],dpi=300)

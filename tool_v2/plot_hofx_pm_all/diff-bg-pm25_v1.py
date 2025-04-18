#!/usr/bin/env python
# plot HofX 
import netCDF4 as nc
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

myvalues=[-20, -10, -5, -3, -1, 1, 3, 5, 10, 20]

f1=nc.Dataset(sys.argv[1],'r')
f2=nc.Dataset(sys.argv[4],'r')
tlon=f1['MetaData/longitude'][:]
tlat=f1['MetaData/latitude'][:]
if ( tlon.max() > 280 ) & (tlon.min() < 180):
 clon=180
else:
 clon=0
 
pm=f1['hofx/particulatematter2p5Surface'][:]-f2['hofx/particulatematter2p5Surface'][:]

if clon==180:
  cproj=ccrs.PlateCarree(central_longitude=180)
  ax=plt.axes(projection=cproj)
  plt.subplots_adjust(wspace=0.6,hspace=0.4)
else:
  ax=plt.axes(projection=ccrs.PlateCarree())
plt.scatter(tlon,tlat,s=5,c=pm,transform=ccrs.PlateCarree(),\
   linewidth=.50,cmap='seismic',norm=colors.BoundaryNorm(boundaries=myvalues,ncolors=256))
plt.axis([-140,-58,24,55])   
cbar=plt.colorbar(fraction=0.07,orientation='horizontal')
cbar.set_label('Surface PM2.5 Diff ($\mu$g/m$^3$)',fontsize=16)
plt.title(sys.argv[2],fontsize=18)
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(states_provinces, edgecolor='gray')
# plt.tight_layout()
plt.savefig(sys.argv[3],dpi=300)

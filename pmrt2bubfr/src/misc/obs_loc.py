import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cartopy.crs as ccrs
import pandas as pd

from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
exit()
ax = plt.axes(projection=ccrs.PlateCarree())   # return GeoAxes object
ax.set_title('Conventional Surface Pressure Location')
ax.coastlines()

#for i in lon:
#   if i>180: i=i-360
#ax.plot(lon,lat,'.b',transform=ccrs.PlateCarree())


# global or not
#ax.set_extent([-130, -60, 20, 55])
#ax.set_global()


# gridlines and labels
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')   # add gridlines to the axes

gl.xlabels_top = False
gl.ylabels_left = False
gl.xlines = True            # draw gridline or not
gl.ylines = True

#gl.xlocator = mticker.FixedLocator([-180, -90, 0, 90, 180 ])
#gl.ylocator = mticker.FixedLocator([-80, -45,0,45, 80 ])

gl.xformatter = LONGITUDE_FORMATTER
gl.yformatter = LATITUDE_FORMATTER

#gl.xlabel_style = {'size': 15, 'color': 'black'}
#gl.xlabel_style = {'color': 'red', 'weight': 'bold'}

# read data
file='2006090812_ges_osse_ps'
data=pd.read_table(file,header=None, delimiter=r"\s+")
lat=data[0]
lon=data[1]
ax.plot(lon,lat,'.b',markersize=1,transform=ccrs.PlateCarree())

# Add labels and title
#ax.set_xlabel('Refractivity (N)')
#ax.set_ylabel('Height (KM)')
#ax.set_title('GPSRO Refractivity')

plt.show()


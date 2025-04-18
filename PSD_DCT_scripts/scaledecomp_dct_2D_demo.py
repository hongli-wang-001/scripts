###################################################################################
# Scale Decomposition
# Based on the python code from Dr. Aaron Johnson (OU/SoM, OU/MAPS)
# (here modified for scale decomposition of analysis increments.)
###################################################################################
#!/bin/env python
import glob, os, sys
import dct_tools
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import matplotlib.axes as maxes
import matplotlib.pyplot as plt
import pprint

################
### Settings ###
################

datadir = './'
filename = 'rtma3d_anl_2019111112.nc'
filename0 = 'rtma3d_fgs_2019111112.nc'

var_name = 'U10'
# var_name = 'V10'
# var_name = 'TH2'
# var_name = 'Q2'
wrffile = datadir + '/' + filename
wrffile0 = datadir + '/' + filename0
domain = 'd01'

# define the scales you want to plot, in km
# what will be plotted is all features > large scale, all features between large and small scale, and all features < small scale
largescale = 400 # km
mesoascale = 100 # km
mesobscale = 40 # km
cnvctscale = 10 # km

L_DETREND = False     # no need to detrend if using DCT

#############################
### Begin executable code ###
#############################

# Routine for plotting WRF data on Basemap
def drawmap(PDATA, TITLE, LOC, outimage):

	F = plt.gcf()
	ax = F.add_subplot(LOC)

        xlon0 = nc.variables['XLONG'][0]
        xlat0 = nc.variables['XLAT'][0]

        print('drawmap: center   lon/lat = ', cen_lon, cen_lat)
        print('drawmap: standard lat1/2  = ',truelat1, truelat2)

        # cropped domain
        print('drawmap: nx, ny (domain cropped) = ',nx,ny)
        xlon  = xlon0[0:ny,0:nx]
        xlat  = xlat0[0:ny,0:nx]

        # using corner-limit mode (so as to match the map to data)
        lon_ll = xlon0[   0,   0]
        lat_ll = xlat0[   0,   0]
        lon_ur = xlon0[ny-1,nx-1]
        lat_ur = xlat0[ny-1,nx-1]

	# Draw the base map behind it with the lats and
	# lons calculated earlier
	m = Basemap(resolution='i',projection='lcc',       \
		width=width_meters,height=height_meters,   \
		lat_0=cen_lat,    lon_0=cen_lon,           \
                lat_1=truelat1,   lat_2=truelat2,          \
		llcrnrlon=lon_ll, llcrnrlat = lat_ll,      \
		urcrnrlon=lon_ur, urcrnrlat = lat_ur)

	# This sets the standard grid point structure at full resolution
	# x,y = m(nc.variables['XLONG'][0],nc.variables['XLAT'][0])
	x,y = m(xlon, xlat)

	m.drawstates(color='k', linewidth=1.25)
	m.drawcoastlines(color='k')
	m.drawcountries(color='k', linewidth=1.25)

        m.drawparallels(np.arange(10,70,20),labels=[1,1,0,0], linewidth=0.4, fontsize=3)
        m.drawmeridians(np.arange(-120,0,20),labels=[0,0,1,1], linewidth=0.4, fontsize=3)

	# PLOT = ax.contourf(x, y, PDATA, LEVS, cmap='seismic',extend='both')
        cm = plt.get_cmap('seismic')
        norm = matplotlib.colors.BoundaryNorm(LEVS, cm.N)
	PLOT = ax.pcolormesh(x, y, PDATA,       cmap=cm, norm=norm               )
	ax.set_title(title, fontsize=8)

	# If final plot, add some other things
	if LOC == 326: 

		#code to make the colorbar outside of the main axis, on the bottom, and lined up
		ax = plt.gcf().add_axes([0.25,-0.03,0.5,0.02])  # Gets the current axes
		bar = plt.colorbar(PLOT,cax=ax,orientation='horizontal',format='%.1f',extend='both') # Plots colorbar in new axis 
		bar.ax.tick_params(labelsize=8)
		bar.ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1.0)) # Make the colorbars numbers nice
		bar.update_ticks()

		plt.tight_layout()
		plt.savefig(outimage,bbox_inches='tight', format='png', dpi=300)


# Set up basemap plotting domain 
nc = Dataset(wrffile, 'r')
nx = len(nc.dimensions['west_east'])
ny = len(nc.dimensions['south_north'])
# nz = len(nc.dimensions['bottom_top'])

dx = float(nc.DX)
dy = float(nc.DY)
width_meters = dx * (nx - 1)
height_meters = dy * (ny - 1)
cen_lat = float(nc.CEN_LAT)
cen_lon = float(nc.CEN_LON)
truelat1 = float(nc.TRUELAT1)
truelat2 = float(nc.TRUELAT2)
standlon = float(nc.STAND_LON)

print('dx = ',dx,' dy =',dy)
print('nx = ',nx,' ny =',ny)

# convert scales to points for 2d DCT
largepoint = largescale / (dx/1000.)
mesoapoint = mesoascale / (dx/1000.)
mesobpoint = mesobscale / (dx/1000.)
cnvctpoint = cnvctscale / (dx/1000.)

# Get the WRF data 
print "Working on %s" % wrffile
time = 0
nc = Dataset(wrffile, 'r')
var_data = nc.variables[var_name][time]
pprint.pprint(var_data.shape)
DATA = var_data[0:ny,0:nx]
print("size of array Data1:",DATA.shape)
if var_name == 'Q2':
  DATA = DATA * 1000.0
print(var_name," of analysis     ----> Max:  %.3f || Min:  %.3f" % (np.ma.max(DATA),np.ma.min(DATA)))

nc0 = Dataset(wrffile0, 'r')
var_data0 = nc0.variables[var_name][time]
pprint.pprint(var_data0.shape)
DATA0 = var_data0[0:ny,0:nx]
if var_name == 'Q2':
  DATA0 = DATA0 * 1000.0
print(var_name," of firstguess   ----> Max:  %.3f || Min:  %.3f" % (np.ma.max(DATA0),np.ma.min(DATA0)))

DATA = DATA - DATA0
print(var_name," of analysis Increments  ----> Max:  %.3f || Min:  %.3f" % (np.ma.max(DATA),np.ma.min(DATA)))

# Get some string info for the time real quick
out = 'ScaleDecomp_of'+'_'+var_name+'_'+'RTMA3D_AnlIncrements'+'.png' # what to call the output images

########################################################
### Do the 2D Discrete Cosine Transform Calculations ###
########################################################

if L_DETREND :
  print('  detrending the raw data ... ')
  SDCT1 = dct_tools.dct2(dct_tools.detrend2d(DATA))
else:
  print('  NO  detrending the raw data ... ')
  SDCT1 = dct_tools.dct2(DATA)

SDCT2 = np.copy(SDCT1)
SDCT3 = np.copy(SDCT1)
SDCT4 = np.copy(SDCT1)
SDCT5 = np.copy(SDCT1)
SDCT6 = np.copy(SDCT1)

# Calculate the 2D wavenumbers at each gridpoint following Denis et al. (2002)
waven = np.zeros_like(DATA)
print("size of array waven:",waven.shape)
for ii in range(0,nx):
	for jj in range(0,ny):
		waven[jj,ii] = ((np.pi*ii/nx)**2 + (np.pi*jj/ny)**2)**0.5

# Filter out certain wavelengths (in wavenumber notation, though) 
for ii in range(0,nx):
	for jj in range(0,ny):

		# Filter out any wavelengths less than the large scale we are interested in 
		# Equvalent to saying: filter out any wavenumbers greater than the large scale we are interested in 
		if waven[jj,ii] > (2*np.pi/largepoint):
			SDCT2[jj,ii] = 0

		# Filter out any wavelengths that are not between the large and small scale
		if waven[jj,ii] > (2*np.pi/mesoapoint) or waven[jj,ii] <= (2*np.pi/largepoint):
			SDCT3[jj,ii] = 0

		# Filter out any wavelengths that are not between the small and mini scale
		if waven[jj,ii] > (2*np.pi/mesobpoint) or waven[jj,ii] <= (2*np.pi/mesoapoint):
			SDCT4[jj,ii] = 0

		# Filter out any wavelengths that are not between the small and mini scale
		if waven[jj,ii] > (2*np.pi/cnvctpoint) or waven[jj,ii] <= (2*np.pi/mesobpoint):
			SDCT5[jj,ii] = 0

		# Filter out any wavelengths greater than the mini scale we are interested in 
		# Equvalent to saying: filter out any wavenumbers less than the mini scale we are interested in 
		if waven[jj,ii] <= (2*np.pi/cnvctpoint):
			SDCT6[jj,ii] = 0


# Rebuild the field using the inverse 2D-DCT
PDATA1 = dct_tools.idct2(SDCT1)
PDATA2 = dct_tools.idct2(SDCT2)
PDATA3 = dct_tools.idct2(SDCT3)
PDATA4 = dct_tools.idct2(SDCT4)
PDATA5 = dct_tools.idct2(SDCT5)
PDATA6 = dct_tools.idct2(SDCT6)
# PDATA7 = PDATA2 + PDATA3 + PDATA4 + PDATA5 + PDATA6

#######################
### Do the Plotting ###
#######################

if var_name == 'U10' or var_name == 'V10':
  LEVS = np.arange(-5,5.1,0.1)
elif var_name == 'TH2':
  LEVS = np.arange(-2,2.1,0.1)
elif var_name == 'Q2':
  LEVS = np.arange(-2,2.1,0.1)

title = 'A) Original function' + ' of ' + var_name + "(Analysis Increments)"
print('----> draw  : '+title)
drawmap(PDATA1, title, 321, out)
	
title = 'B) Wavelengths >= %i km only' % largescale
print('----> draw  : '+title)
drawmap(PDATA2, title, 322, out)

title = 'C) Wavelengths >= %i km and < %i km only' % (mesoascale, largescale)
print('----> draw  : '+title)
drawmap(PDATA3, title, 323, out)

title = 'D) Wavelengths >= %i km and < %i km only' % (mesobscale, mesoascale)
print('----> draw  : '+title)
drawmap(PDATA4, title, 324, out)

title = 'E) Wavelengths >= %i km and < %i km only' % (cnvctscale, mesobscale)
print('----> draw  : '+title)
drawmap(PDATA5, title, 325, out)

title = 'F) Wavelengths < %i km only' % cnvctscale
print('----> draw  : '+title)
drawmap(PDATA6, title, 326, out)

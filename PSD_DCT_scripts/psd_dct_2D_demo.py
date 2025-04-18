#################################################################################
#        Compute power spectral density (PSD) with                              #
#         Discrete Cosine/Sine Transform (DCST)                                 #
#                                                                               #
# NOAA/NWS/NCEP/EMCr,RTMA Group, Gang Zhao                                      #
# Code History:                                                                 #
#              Originally adopting some code from Dr. Jacob Carley (NCEP/EMC)   #
#              and Dr.Youngsun Jung (OU/CAPS)                                   #        
# Question: email Gang.Zhao@NOAA.GOV                                            #
#################################################################################
#!/bin/env python

from __future__ import print_function
import sys,os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from numpy import array
# import scipy.fftpack
# import time as mytime
# from datetime import datetime,timedelta
from netCDF4 import Dataset
import pygrib

# sys.path.append('/home/Gang.Zhao/local/python/dct_tools')
# import dct_tools
import dcst

######################################################################
def wn2wl(x):
  return np.pi*2./(x*1000.)

######################################################################
def wl2wn(x):
  return np.pi*2./(x*1000.)

######################################################################
def ke_2d(u):
  nshp = np.shape(u)
  print('ke_2d: shape of u -->',nshp) 
  ny = nshp[0]
  nx = nshp[1]
  print('ke_2d: nx ny -->',nx,ny) 
  ke = 0.0
  for i in np.arange(nx):
    for j in np.arange(ny):
      ke = ke + u[j,i]*u[j,i]
  return ke

######################################################################

def dct_psd_2d(u,nx_in,ny_in,dx_in,dy_in):

  xi = 0
  yi = 0
  nx = nx_in
  ny = ny_in
  domainx=nx                # verification domain
  domainy=ny                # verification domain
  xf = xi+domainx
  yf = yi+domainy
  dx, dy = (dx_in,dy_in)    # Unit: m
  grid_spacing = (dy, dx)
  nx,ny = (domainx,domainy)
  dims = (ny,nx)
  scale_factor = 1.e3    # For m -> km
  print("xi yi:",xi,yi,"  xf yf:",xf,yf,"  nx ny:",nx,ny)

#----------------------------#
# Initialization of arrays   #
#----------------------------#
  utemp_var = np.zeros((ny,nx))
  vtemp_var = np.zeros((ny,nx))
  total_Kkv = np.zeros((ny,nx))

# wavenumbers
  L = (ny*dy, nx*dx)
  print("domain size: L are " ,  L, ' m')
  kx0 = 2.*np.pi*np.fft.fftfreq(2*nx)*nx/L[1]
  kx  = kx0[0:nx]
  ky0 = 2.*np.pi*np.fft.fftfreq(2*ny)*ny/L[0]
  ky  = ky0[0:ny]

  k_x = kx*L[1]
  k_y = ky*L[0]

# physical limits to the wavenumbers
  kmin = 2.*np.pi*1.0/np.min([L[0],L[1]])
  kmax = 2.*np.pi*np.min([0.5*dims[0]/L[0],0.5*dims[1]/L[1]])
  print("kmin, kmax:",kmin, kmax)

  kbins = np.arange(kmin, kmax, kmin)
  N = len(kbins)

# bin the Fourier KE into radial kbins
  ky2d, kx2d = np.meshgrid(ky, kx, indexing="ij")
  k = np.sqrt(kx2d**2 + ky2d**2 )
# print("shape of array k:")
# print(np.shape(k))

  whichbin = np.digitize(k.flat, kbins)
# print("shape of array whichbin:")
# print(np.shape(whichbin))
  
  ncount = np.bincount(whichbin)

  k = 0.5*(kbins[0:N-1] + kbins[1:N])
# print("shape of array k (after kbins):")
# print(np.shape(k))

  wavelen = 2.*np.pi/k/scale_factor   # m to km 

  total_power = np.zeros(len(ncount)-2)

  utemp_var = u[yi:yf,xi:xf]
  print("shape of array of variable (u):",u.shape,"  shape of utemp_var:",utemp_var.shape)

# Kenetic Energy in physical space
  ke_u = ke_2d(utemp_var)
  print('dct_psd_2d:  ke_u (from ke_2d) -->',ke_u)
  ke_total = ke_u

# FFT/DCST
  Kkkv_u = dcst.dct2(utemp_var)
  temp_kvfft = np.abs(Kkkv_u*np.conj(Kkkv_u)                         )
  temp_kvfft = temp_kvfft / (2.*nx*2.*ny) # scaling (no averaging)
  temp_kvfft[0,:] = temp_kvfft[0,:] / 2.  # further scaling for 0-wave in one-direction
  temp_kvfft[:,0] = temp_kvfft[:,0] / 2.  # further scaling for 0-wave in another-direction
  ke_from_fft = np.real(np.sum(temp_kvfft))
  ke_ratio=ke_from_fft/ke_total
  print(" ----> checking Parceval Identity: " )
  print(" --------> KE= %.8f  kefft=%.8f  kefft/KE= %.8f" % (ke_total,ke_from_fft,ke_ratio))

# scaled/normalized, domain-averaged, power-spectrum-density (PSD)
  total_Kkv = 0.5*np.abs(Kkkv_u*np.conj(Kkkv_u))/((2*nx*2*ny)*(nx*ny)*kmin)
  total_Kkv[0,:] = total_Kkv[0,:] / 2.
  total_Kkv[:,0] = total_Kkv[:,0] / 2.

# Energy Spectra (from 2-D wavenumber to 1-D wavenumber)
  E_spectrum = np.zeros(len(ncount)-1)

  for n in range(1,len(ncount)):
      E_spectrum[n-1] = np.sum(total_Kkv.flat[whichbin==n])

  for n in range (1, len(ncount)-1):
      total_power[n-1] = E_spectrum[n-1]
#     print(num_states, n, power(n), total_power(n))

# Compute max energy
  index = np.argmax(E_spectrum)
  kmax = k[index]
  Emax = E_spectrum[index]
  print("index = ", index, "  Emax = ", Emax, "  kmax = ", kmax)
  print("size of spectrum:", len(E_spectrum))

  return total_power, k

######################################################################
# Main Driver Program to test the function of dct_psd_2D

if __name__ == '__main__':

  print('Starting...')

  dx = 3000.
  dy = 3000.
  DX = dx

#=================================================================================#
  adate = '2019-11-11_12_00_00 UTC'
  var_name = 'TH2'                      # U10; V10; TH2; Q2; # WIND10(not working now);
  subset = var_name
  print("subset:",subset)
  if var_name == 'TH2':
      vname_grb = '2 metre temperature'
      vname     = var_name
  elif var_name == 'Q2':
      vname_grb = 'Specific humidity'
      vname     = var_name
  elif var_name == 'U10':
      vname_grb = '10 metre U wind component'
      vname     = var_name
  elif var_name == 'V10':
      vname_grb = '10 metre V wind component'
      vname     = var_name
  elif var_name == 'WIND10':
      vname_grb = '10 metre wind speed'
  datahome = './'

#=================================================================================#
# get 2D data -- rtma3D in netcdf format
#=================================================================================#
  datadir= datahome
  filename='rtma3d_anl_2019111112.nc'
  print('  ---> processing data file (netcdf RTMA-3D):',datadir,'/',filename,sep='')
# if grib
#  filegrp_grbs2 = pygrib.open(datadir+'/'+filename)        # pygrib
#  print('----> information about the analysis grib2 file: ')
#  print('    file name : ----->',filegrp_grbs2.name)
#  print('    number of total messages in file : ----->',filegrp_grbs2.messages)
#  filegrp_grbs2.seek(0)
#  for grb in filegrp_grbs2 :
#      print('    mesg of grb:',grb)
#  grbmsg = filegrp_grbs2.select(name=vname_grb)[0]
#  DATA2D = grbmsg.values
#  print('    ----> data2D array size:',DATA2D.shape)
#  print('    ----> data2D max and min:',DATA2D.max(), DATA2D.min())
# if netcdf
  filegrp = Dataset(datadir+'/'+filename, 'r')             # netCDF4
  DATA2D = filegrp.variables[vname][0,:,:]
  if var_name == 'Q2':
      DATA2D = DATA2D * 1000.0           # kg/kg --> g/kg
  print('    ----> data2D(from netcdf) array size:',DATA2D.shape)
  print("Variable   ",vname," --> Max :  %.3f ; Min :  %.3f" % (np.ma.max(DATA2D),np.ma.min(DATA2D)))
#=================================================================================#
# Power Spectram Analysis for variables
  ny_lg,nx_lg = np.shape(DATA2D)
  print('dimension size: ',nx_lg, ny_lg)
  # FFT
  psd,wn = dct_psd_2d(DATA2D,nx_lg,ny_lg,dx,dy)
  wl = 2.0 * np.pi / (wn*1000.)

#=================================================================================#
# plot
#=================================================================================#
  fig   = plt.figure(figsize=(6,6), facecolor='w', dpi=300)
  ax1	= fig.add_subplot(111)             # python 2.x
# ax1	= fig.subplots()                   # python 3.x

  anl, = ax1.plot(wn, psd, color='r',lw=2)

  lindborg = (9.1e-4*np.power(wn, -5./3.) + 3e-10*np.power(wn, -3.)) # this is correct based on Fig. 7 in Lindborg et al. (1999); not sure how valid below tropopause
  ymin = psd.min()
  ymax = psd.max()
  lind, = ax1.plot(wn, lindborg, color='k', lw=1, linestyle='-.')

# ax1.grid()
# ax1.grid(which='major', axis='y')
  ax1.set_axisbelow(True)
  ax1.loglog()

  if var_name == 'U10' or var_name == 'V10' or var_name == 'WIND10':
      ax1.set_ylabel(r'PSD:  E($\kappa$)  [m$^3$ s$^{-1}$]')     # U10 V10 WIND10
  if var_name == 'TH2':
      ax1.set_ylabel(r'PSD:  E($\kappa$)  [K$^2$ m$^1$]')        # TH2
  if var_name == 'Q2':
      ax1.set_ylabel(r'PSD:  E($\kappa$)  [(g/kg)$^2$ m$^1$]')    # Q2

  ax1.set_xlabel(r'Wavenumber [radians m$^{-1}$]')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# adding secondary x-axis for wavelength
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   the following code only works with matplotlib 3.1.0 and beyond.
# ax2 = ax1.secondary_xaxis('top',functions = (wn2wl, wl2wn))
# ax2.set_xlabel(b'Wavelength[km]')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
#   alternative way to add 2nd axis to the plot if no matplotlib 3.1.0 and beyond
  ax2 = ax1.twiny()
# if loglog used for ax1 and ax2, the following two lines must be done before set _xticks
  ax2.set_xlim(ax1.get_xlim())			# x-axis 2 matches x-axis 1, so the xticks could correspondingly be located.
# ax2.set_xbound(ax1.get_xbound())
  ax2.loglog()					# twiny does not make ax2 automatically following the log-scale of ax1.
						# if using secondary_xaxis, the log-scale could be inherited to children class.
# ax2.grid()
# ax2.grid(which='major', axis='x')
  ax2.set_axisbelow(True)

  ax2_labels	= [5., 10., 50., 100., 500., 1000.]	# labels for 2nd x-axis (wavelength in km)
  wl_to_wn	= lambda wl: 2.0*np.pi/(wl*1000.0)
  ax2_lb_pos	= [wl_to_wn(wl) for wl in ax2_labels]
  ax2.set_xticks(ax2_lb_pos)
  ax2.set_xticklabels(ax2_labels, rotation='45',fontsize='xx-small')

# ax2.xaxis.set_ticks_position('bottom')
# ax2.xaxis.set_label_position('bottom')
# ax2.spines['bottom'].set_position(('outward',36))

  ax2.set_xlabel(r'Wavelength [km]', fontsize='xx-small')

  ax1.legend([anl, lind],['analysis_rtma3D', 'lindborg:k$^{-5/3}$ + k$^{-3}$ '], fontsize='xx-small')

  fig_title = ax1.set_title('Power Spectrum Density of Analysis of RTMA-3D (HRRR domain) '+' \n '+var_name+'   '+adate+' (using DCT)')
  fig_title.set_y(1.15)

  fig.subplots_adjust(top=0.85)

  fig.savefig('PSD_of_'+var_name+'.png', bbox_inches='tight', dpi=600) 

  plt.clf()
  plt.close()

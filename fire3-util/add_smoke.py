import os
import netCDF4 as nc
from netCDF4 import Dataset

#Add smoke and dust to re-start file

file_to_extract= "tracer.res.tile1.nc"
open_fext=nc.Dataset(file_to_extract)
smoke_2_add=open_fext['smoke'][0,:,:,:]
print('SMOKE MAX to ADD',smoke_2_add.max())
dust_2_add=open_fext['dust'][0,:,:,:]
dust_3_add=open_fext['coarsepm'][0,:,:,:]
file='gfs_data.tile7.halo0.nc'
file_input= nc.Dataset(file,'r+',format='NETCED4')
smoke=file_input.createVariable('smoke', 'f8',('lev','lat', 'lon'))
smoke[:,:,:]=0
smoke[1:66,:,:]=smoke_2_add
smoke.units = " ug/kg"
dust=file_input.createVariable('dust', 'f8',('lev','lat', 'lon'))
dust[:,:,:]=0
dust[1:66,:,:]=dust_2_add
dust.units = " ug/kg"
coarsepm=file_input.createVariable('coarsepm', 'f8',('lev','lat', 'lon'))
coarsepm[:,:,:]=0
coarsepm[1:66,:,:]=dust_3_add
coarsepm.units = " ug/kg"
print('SMOKE MAX ADDED',smoke.max)
print('SMOKE MAX ADDED',dust.max)
print('SMOKE MAX ADDED',coarsepm.max)
file_input.variables
file_input.close()



#######################################################################

import os, glob
import numpy as np
import datetime as dt
import netCDF4 as nc
import time
import pickle


def addmodelvar_derived(restartpath):

  ########################################################################
  #                                                                      #
  # add "var_derived" variable to control member and ensemble restart phys files #
  #                                                                      #
  ########################################################################

  # specify some variables used to convert pm2p5
  scale_ug2kg=1.0e-9

  # make variables for names of model files
  corefile = restartpath+'fv_core.res.tile1.nc' # need diff names for ensembles?
  tracerfile = restartpath+'fv_tracer.res.tile1.nc'
  physfile = restartpath+'tmp.nc'

  # open core file to read in delta pressure, then close
  u = nc.Dataset(corefile,'r')
  delp = u['delp'][0] # delta pressure in Pascals
  u.close()

  # open tracer file to read pm2p5, then close
  u = nc.Dataset(tracerfile,'r')
  pm2p5 = u['pm25_tot'][0] # units kg/kg 
  u.close()

  # open physics file to write column pm2.5 total 
  u = nc.Dataset(physfile,'r+')

  #if 'column_pm2p5_total' in u.variables.keys():
  #  print('column_pm2p5_total field exists! Exiting')
  #  u.close()
  #  return()


  # units of [delp/g] = Pa*s^2/m = N*s^2/m^3 = kg/m^2 = dry air mass per m^2 in the cell
  # * by 3000 m * 3000 m = 9e6=.009e9 for air mass in cell, in kg
  # then * by pm2p5 in kg/kg to get total pm2p5 mass in cell
  # and sum over the column for pm2p5 mass in column kg/m^2*kg/kg = kg/m^2
  print (delp.shape)
  print (pm2p5.shape)
  pm2p5_per_area = np.sum(scale_ug2kg*delp*pm2p5/9.8,axis=0)

  var_derived = pm2p5_per_area

  # add z- and Time dimensions
  var_derived = np.array([[var_derived for i in range(u.dimensions['zaxis_1'].size)]],np.float32)

  # add 4D column_pm2p5_total grid to tracer NetCDF and close
  # Time,zaxis_1,yaxis_1,xaxis_1 are dimensions to write
  if 'column_pm2p5_total' in u.variables.keys():
    print('Overwriting column_pm2p5_total values.')
    u['column_pm2p5_total'][:] = var_derived.astype('float32') 
  else:
    var_derived_out = u.createVariable('column_pm2p5_total',np.float32,('Time','zaxis_1','yaxis_1','xaxis_1'),chunksizes=(1,1,u.dimensions['yaxis_1'].size,u.dimensions['xaxis_1'].size))
    var_derived_out[:] = var_derived.astype('float32')
  u.close()
  return()

if __name__=="__main__":

  ########################################################
  #                                                      #
  # add model column_pm2p5_total to input or ensemble files in place    #
  #                                                      #
  ########################################################
    nwges_dir = os.environ.get("nwges_dir")
    restart_prefix = os.environ.get("restart_prefix")
    this_path = nwges_dir+'/RESTART/'+restart_prefix+'.'
    print('Model data at '+this_path+'*')
#    if os.path.exists(this_path+'tmp.nc'):
    addmodelvar_derived(this_path)
#    else:
#      print('Model data not found at '+this_path+'*')

    quit()


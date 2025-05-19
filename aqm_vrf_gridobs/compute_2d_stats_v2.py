import numpy as np
import netCDF4 as nc
import glob
import os

# Find all relevant NetCDF files
#file_list = sorted(glob.glob('vrf_3hcyc_pm_h100v12/*00/pm25.*_f001.nc'))
file_list = sorted(glob.glob('vrf_3hcyc_ctr/*00/pm25.*_f001.nc'))

if not file_list:
    raise FileNotFoundError("No matching pm25.*_f001.nc files found.")

# Stacks for pm25_diff and pm25
pm25_diff_stack = []
pm25_stack = []
pm25_aqm_stack = []

for file in file_list:
    with nc.Dataset(file, 'r') as ds:
        # Load pm25_diff
        if 'pm25_diff' not in ds.variables:
            raise KeyError(f"'pm25_diff' not found in {file}")
        pm25_diff = ds.variables['pm25_diff'][:]
        pm25_diff = np.where(pm25_diff == -1, np.nan, pm25_diff)
        pm25_diff_stack.append(pm25_diff)

        # Load pm25 obs 
        if 'pm25' not in ds.variables:
            raise KeyError(f"'pm25' not found in {file}")
        pm25 = ds.variables['pm25'][:]
        pm25 = np.where(pm25 == -1, np.nan, pm25)
        pm25_stack.append(pm25)

        # Load pm25_aqm  model prediction 
        if 'pm25_aqm' not in ds.variables:
            raise KeyError(f"'pm25_aqm' not found in {file}")
        pm25_aqm = ds.variables['pm25_aqm'][:]
        pm25_aqm = np.where(pm25_aqm == -1, np.nan, pm25_aqm)
        pm25_aqm_stack.append(pm25_aqm)

# Convert to 3D arrays [file, y, x]
pm25_diff_array = np.array(pm25_diff_stack)
pm25_array = np.array(pm25_stack)
pm25_aqm_array = np.array(pm25_aqm_stack)

# Compute stats
pm25_diff_mean = np.nanmean(pm25_diff_array, axis=0)
pm25_diff_rms = np.sqrt(np.nanmean(pm25_diff_array**2, axis=0))
pm25_mean = np.nanmean(pm25_array, axis=0)
pm25_aqm_mean = np.nanmean(pm25_aqm_array, axis=0)

# Get dimensions from first file
with nc.Dataset(file_list[0], 'r') as ds_in:
    y_dim = ds_in.dimensions['y'].size
    x_dim = ds_in.dimensions['x'].size

# Create output NetCDF
with nc.Dataset('pm25_vrfgrid_stats.nc', 'w', format='NETCDF4') as ds_out:
    ds_out.createDimension('y', y_dim)
    ds_out.createDimension('x', x_dim)

    # Variables
    diff_mean_var = ds_out.createVariable('pm25_diff_mean', 'f4', ('y', 'x'), zlib=True)
    diff_rms_var = ds_out.createVariable('pm25_diff_rms', 'f4', ('y', 'x'), zlib=True)
    pm25_mean_var = ds_out.createVariable('pm25_obs__mean', 'f4', ('y', 'x'), zlib=True)
    pm25_aqm_mean_var = ds_out.createVariable('pm25_aqm_mean', 'f4', ('y', 'x'), zlib=True)

    # Assign data
    diff_mean_var[:, :] = pm25_diff_mean
    diff_rms_var[:, :] = pm25_diff_rms
    pm25_mean_var[:, :] = pm25_mean
    pm25_aqm_mean_var[:, :] = pm25_aqm_mean

    # Add metadata
    diff_mean_var.units = 'µg/m³'
    diff_mean_var.description = 'Mean of PM2.5 differences (pm25_diff)'
    
    diff_rms_var.units = 'µg/m³'
    diff_rms_var.description = 'Root Mean Square of PM2.5 differences (pm25_diff)'

    pm25_mean_var.units = 'µg/m³'
    pm25_mean_var.description = 'Mean of PM2.5 OBS values'

    pm25_aqm_mean_var.units = 'µg/m³'
    pm25_aqm_mean_var.description = 'Mean of PM2.5 FCST values'
print("Saved statistics to pm25_vrfgrid_stats.nc")

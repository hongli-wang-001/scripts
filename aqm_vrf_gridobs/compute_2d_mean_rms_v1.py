import numpy as np
import netCDF4 as nc
import glob
import os

# Find all relevant NetCDF files
file_list = sorted(glob.glob('vrf_3hcyc_pm_h100v12/*00/pm25.*_f001.nc'))

if not file_list:
    raise FileNotFoundError("No matching pm25.*_f001.nc files found.")

# Read data and stack
pm25_diff_stack = []

for file in file_list:
    with nc.Dataset(file, 'r') as ds:
        if 'pm25_diff' not in ds.variables:
            raise KeyError(f"'pm25_diff' not found in {file}")
        pm25_diff = ds.variables['pm25_diff'][:]
        pm25_diff = np.where(pm25_diff == -1, np.nan, pm25_diff)  # mask missing values
        pm25_diff_stack.append(pm25_diff)

# Convert to 3D array [file, y, x]
pm25_diff_array = np.array(pm25_diff_stack)

# Compute mean and RMS across files, ignoring NaNs
pm25_mean = np.nanmean(pm25_diff_array, axis=0)
pm25_rms = np.sqrt(np.nanmean(pm25_diff_array**2, axis=0))

# Use shape info from first file
with nc.Dataset(file_list[0], 'r') as ds_in:
    y_dim = ds_in.dimensions['y'].size
    x_dim = ds_in.dimensions['x'].size

# Create output NetCDF
with nc.Dataset('pm25_diff_2d_stats.nc', 'w', format='NETCDF4') as ds_out:
    # Dimensions
    ds_out.createDimension('y', y_dim)
    ds_out.createDimension('x', x_dim)

    # Variables
    mean_var = ds_out.createVariable('pm25_diff_mean', 'f4', ('y', 'x'), zlib=True)
    rms_var = ds_out.createVariable('pm25_diff_rms', 'f4', ('y', 'x'), zlib=True)

    # Assign data
    mean_var[:, :] = pm25_mean
    rms_var[:, :] = pm25_rms

    # Attributes
    mean_var.units = 'µg/m³'
    rms_var.units = 'µg/m³'
    mean_var.description = 'Mean of PM2.5 differences'
    rms_var.description = 'Root Mean Square of PM2.5 differences'

print(f"Saved statistics to pm25_diff_2d_stats.nc")

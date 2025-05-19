import numpy as np
import netCDF4 as nc
import glob
import os

for hour in range(1, 25):
    fhr = f"f{hour:03d}"
    file_list = sorted(glob.glob(f"vrf_3hcyc_pm_h100v12/*00/pm25.*_{fhr}.nc"))
    #file_list = sorted(glob.glob(f"vrf_3hcyc_ctr/*00/pm25.*_{fhr}.nc"))

    if not file_list:
        print(f"No files found for {fhr}, skipping.")
        continue

    print(f"Processing {len(file_list)} files for {fhr}...")

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

    pm25_diff_array = np.array(pm25_diff_stack)
    pm25_array = np.array(pm25_stack)
    pm25_aqm_array = np.array(pm25_aqm_stack)

    pm25_diff_mean = np.nanmean(pm25_diff_array, axis=0)
    pm25_diff_rms = np.sqrt(np.nanmean(pm25_diff_array**2, axis=0))
    pm25_mean = np.nanmean(pm25_array, axis=0)
    pm25_aqm_mean = np.nanmean(pm25_aqm_array, axis=0)

    # Get dimension sizes from the first file
    with nc.Dataset(file_list[0], 'r') as ds_in:
        y_dim = ds_in.dimensions['y'].size
        x_dim = ds_in.dimensions['x'].size

    # Write stats to NetCDF
    out_file = f"pm25_2d_stats_{fhr}.nc"
    with nc.Dataset(out_file, 'w', format='NETCDF4') as ds_out:
        ds_out.createDimension('y', y_dim)
        ds_out.createDimension('x', x_dim)

        mean_var = ds_out.createVariable('pm25_diff_mean', 'f4', ('y', 'x'), zlib=True)
        rms_var = ds_out.createVariable('pm25_diff_rms', 'f4', ('y', 'x'), zlib=True)
        pm25_mean_var = ds_out.createVariable('pm25_obs__mean', 'f4', ('y', 'x'), zlib=True)
        pm25_aqm_mean_var = ds_out.createVariable('pm25_aqm_mean', 'f4', ('y', 'x'), zlib=True)

        mean_var[:, :] = pm25_diff_mean
        rms_var[:, :] = pm25_diff_rms
        pm25_mean_var[:, :] = pm25_mean
        pm25_aqm_mean_var[:, :] = pm25_aqm_mean

        mean_var.units = 'µg/m³'
        mean_var.description = f'Mean of PM2.5 differences ({fhr})'

        rms_var.units = 'µg/m³'
        rms_var.description = f'RMS of PM2.5 differences ({fhr})'

        pm25_mean_var.units = 'µg/m³'
        pm25_mean_var.description = f'Mean PM2.5 OBS values ({fhr})'

        pm25_aqm_mean_var.units = 'µg/m³'
        pm25_aqm_mean_var.description = 'Mean of PM2.5 FCST values ({fhr})'
    print(f"Saved statistics to {out_file}")

import numpy as np
import netCDF4 as nc
import glob
import os

for hour in range(0, 24):
    fhr = f"f{hour:03d}"
    #file_list = sorted(glob.glob(f"vrf_3hcyc_pm_h100v12/*??/aod.*_{fhr}.nc"))
    #file_list = sorted(glob.glob(f"vrf_3hcyc_ctr/*18/aod.*_{fhr}.nc"))
    file_list = sorted(glob.glob(f"EXP/*HH/aod.*_{fhr}.nc"))
    if not file_list:
        print(f"No files found for {fhr}, skipping.")
        continue

    print(f"Processing {len(file_list)} files for {fhr}...")

    aod_diff_stack = []
    aod_stack = []
    aod_aqm_stack = []

    for file in file_list:
        with nc.Dataset(file, 'r') as ds:
            # Load aod_diff
            if 'aod_diff' not in ds.variables:
                raise KeyError(f"'aod_diff' not found in {file}")
            aod_diff = ds.variables['aod_diff'][:]
            aod_diff = np.where(aod_diff == -1, np.nan, aod_diff)
            aod_diff_stack.append(aod_diff)

           # Load aod obs 
            if 'aod550nm' not in ds.variables:
                raise KeyError(f"'aod' not found in {file}")
            aod = ds.variables['aod550nm'][:]
            aod = np.where(aod == -1, np.nan, aod)
            aod_stack.append(aod)

            # Load aod_aqm  model prediction 
            if 'aod_aqm' not in ds.variables:
                raise KeyError(f"'aod_aqm' not found in {file}")
            aod_aqm = ds.variables['aod_aqm'][:]
            aod_aqm = np.where(aod_aqm == -1, np.nan, aod_aqm)
            aod_aqm_stack.append(aod_aqm)

    aod_diff_array = np.array(aod_diff_stack)
    aod_array = np.array(aod_stack)
    aod_aqm_array = np.array(aod_aqm_stack)

    aod_diff_mean = np.nanmean(aod_diff_array, axis=0)
    aod_diff_rms = np.sqrt(np.nanmean(aod_diff_array**2, axis=0))
    aod_mean = np.nanmean(aod_array, axis=0)
    aod_aqm_mean = np.nanmean(aod_aqm_array, axis=0)

    # Get dimension sizes from the first file
    with nc.Dataset(file_list[0], 'r') as ds_in:
        y_dim = ds_in.dimensions['y'].size
        x_dim = ds_in.dimensions['x'].size

    # Write stats to NetCDF
    out_file = f"aod_2d_stats_{fhr}.nc"
    with nc.Dataset(out_file, 'w', format='NETCDF4') as ds_out:
        ds_out.createDimension('y', y_dim)
        ds_out.createDimension('x', x_dim)

        mean_var = ds_out.createVariable('aod_diff_mean', 'f4', ('y', 'x'), zlib=True)
        rms_var = ds_out.createVariable('aod_diff_rms', 'f4', ('y', 'x'), zlib=True)
        aod_mean_var = ds_out.createVariable('aod_obs_mean', 'f4', ('y', 'x'), zlib=True)
        aod_aqm_mean_var = ds_out.createVariable('aod_aqm_mean', 'f4', ('y', 'x'), zlib=True)

        mean_var[:, :] = aod_diff_mean
        rms_var[:, :] = aod_diff_rms
        aod_mean_var[:, :] = aod_mean
        aod_aqm_mean_var[:, :] = aod_aqm_mean

        mean_var.units = ''
        mean_var.description = f'Mean of AOD differences ({fhr})'

        rms_var.units = ''
        rms_var.description = f'RMS of AOD differences ({fhr})'

        aod_mean_var.units = ''
        aod_mean_var.description = f'Mean AOD OBS values ({fhr})'

        aod_aqm_mean_var.units = ''
        aod_aqm_mean_var.description = 'Mean of AOD FCST values ({fhr})'
    print(f"Saved statistics to {out_file}")

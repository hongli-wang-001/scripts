import netCDF4 as nc
import numpy as np

# Open fcst.nc (read-only) and pm25.nc (append mode)
with nc.Dataset('fcst.nc', 'r') as fcst_ds, nc.Dataset('pm25.nc', 'a') as pm25_ds:
    # Read the specific slice from fcst.nc
    pm25_ave = fcst_ds.variables['pm25_ave'][0, 63, :, :]  # Shape: [225, 393]

    # Read pm25 from pm25.nc
    pm25 = pm25_ds.variables['pm25'][:, :]  # Shape: [225, 393]

    # Compute the difference
    diff = pm25_ave - pm25

    # Remove or mask invalid values if needed (e.g., NaNs)
    valid_diff = diff[~np.isnan(diff)]
 
    # Calculate statistics
    mean_diff = np.mean(valid_diff)
    median_diff = np.median(valid_diff)
    std_diff = np.std(valid_diff)
    rms_diff = np.sqrt(np.mean(valid_diff ** 2))

    # Print statistics
    print("Statistics for pm25_diff (pm25_ave[0,63,:,:] - pm25):")
    print(f"  Mean      : {mean_diff:.4f}")
    print(f"  Median    : {median_diff:.4f}")
    print(f"  Std Dev   : {std_diff:.4f}")
    print(f"  RMS       : {rms_diff:.4f}")

    # If variable already exists, delete it to avoid conflict
    if 'pm25_diff' in pm25_ds.variables:
        del pm25_ds.variables['pm25_diff']

    # Define the new variable in pm25.nc
    diff_var = pm25_ds.createVariable('pm25_diff', 'f4', ('y', 'x'))
    diff_var[:, :] = diff

    # Add metadata
    diff_var.units = 'µg/m³'
    diff_var.long_name = 'Difference: pm25_ave[0,63,:,:] - pm25'

print("Difference written to pm25.nc as variable 'pm25_diff'.")

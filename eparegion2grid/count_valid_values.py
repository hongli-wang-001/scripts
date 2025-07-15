import netCDF4 as nc
import numpy as np

def count_non_nan_values(nc_file_path, var_name):
    # Open the NetCDF file
    with nc.Dataset(nc_file_path, 'r') as ds:
        # Read the variable data
        data = ds.variables[var_name][:]

        # Count non-NaN values
        non_nan_count = np.count_nonzero(~np.isnan(data))

        print(f"Number of non-NaN values in variable '{var_name}': {non_nan_count}")

# Example usage
count_non_nan_values('pm25.nc', 'epa_region')


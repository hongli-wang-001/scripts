import numpy as np
from netCDF4 import Dataset
import glob
import re

# Define file pattern and exclude pattern
file_pattern = 'con/aqm.202009*/18/aqm.t18z.dyn.f001.nc'
exclude_pattern = re.compile(r'f000')

# List files matching the pattern
file_list = sorted(glob.glob(file_pattern))

# Exclude files with the f000 pattern
file_list = [file for file in file_list if not exclude_pattern.search(file)]

# Check if any files are left after exclusion
if not file_list:
    raise ValueError("No files found after excluding the specified pattern.")


# Load dimensions from the first file
first_file = file_list[0]
with Dataset(first_file, 'r') as nc:
    # Get dimensions and shape from pm25_tot
    pm25_tot_dims = nc.variables['pm25_tot'].dimensions
    dimensions = {dim: len(nc.dimensions[dim]) for dim in pm25_tot_dims}
    dim_names = list(dimensions.keys())

    # Print each dimension's name and size
    print("Dimensions and sizes:")
    for dim_name in dim_names:
        dim_size = dimensions[dim_name]
        print(f"{dim_name}: {dim_size}")

# Initialize lists to collect data
data_pm25_tot = []
data_pm25_ave = []

# Load data from each file
for file in file_list:
    with Dataset(file, 'r') as nc:
        if 'pm25_tot' in nc.variables and 'pm25_ave' in nc.variables:
            data_pm25_tot.append(nc.variables['pm25_tot'][:])
            data_pm25_ave.append(nc.variables['pm25_ave'][:])
        else:
            print(f"Missing pm25_tot or pm25_ave in file {file}")

# Convert lists to numpy arrays for easier manipulation
data_pm25_tot = np.concatenate(data_pm25_tot, axis=0)  # Shape: (n_files, ... ) depending on dimensions
data_pm25_ave = np.concatenate(data_pm25_ave, axis=0)  # Shape: (n_files, ... ) depending on dimensions

# Calculate statistics
mean_pm25_tot = np.mean(data_pm25_tot, axis=0)  # Shape: (dimensions excluding time)
std_pm25_tot = np.std(data_pm25_tot, axis=0)
mean_pm25_ave = np.mean(data_pm25_ave, axis=0)
std_pm25_ave = np.std(data_pm25_ave, axis=0)

# Calculate the mean and std of the difference between pm25_tot and pm25_ave
difference = data_pm25_tot - data_pm25_ave
mean_diff = np.mean(difference, axis=0)
std_diff = np.std(difference, axis=0)

# Write results to a new NetCDF file
output_file = 'output_statistics.nc'
with Dataset(output_file, 'w', format='NETCDF4') as nc_out:
    # Create dimensions
    for dim_name in dim_names:
        nc_out.createDimension(dim_name, dimensions[dim_name])

    # Create variables
    nc_out.createVariable('mean_pm25_tot', np.float32, dim_names)
    nc_out.createVariable('std_pm25_tot', np.float32, dim_names)
    nc_out.createVariable('mean_pm25_ave', np.float32, dim_names)
    nc_out.createVariable('std_pm25_ave', np.float32, dim_names)
    nc_out.createVariable('mean_diff', np.float32, dim_names)
    nc_out.createVariable('std_diff', np.float32, dim_names)

    # Write data to variables
    nc_out.variables['mean_pm25_tot'][:] = mean_pm25_tot
    nc_out.variables['std_pm25_tot'][:] = std_pm25_tot
    nc_out.variables['mean_pm25_ave'][:] = mean_pm25_ave
    nc_out.variables['std_pm25_ave'][:] = std_pm25_ave
    nc_out.variables['mean_diff'][:] = mean_diff
    nc_out.variables['std_diff'][:] = std_diff

print("Statistics written to output_statistics.nc")


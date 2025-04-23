import numpy as np
from netCDF4 import Dataset
import glob
import re

# File2-File1; 04/16/2025

# Define file pattern and exclude pattern
file_pattern = 'con/aqm.202009*/00/aqm.t00z.dyn.f001.nc'
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
    # Get dimensions and shape from pm25_ave
    pm25_tot_dims = nc.variables['pm25_ave'].dimensions
    dimensions = {dim: len(nc.dimensions[dim]) for dim in pm25_tot_dims}
    dim_names = list(dimensions.keys())

    # Print each dimension's name and size
    print("Dimensions and sizes:")
    for dim_name in dim_names:
        dim_size = dimensions[dim_name]
        print(f"{dim_name}: {dim_size}")

# Define file pattern and exclude pattern
file_pattern2 = 'pmf/aqm.202009*/00/aqm.t00z.dyn.f001.nc'

# List files matching the pattern
file_list2 = sorted(glob.glob(file_pattern2))

# Exclude files with the f000 pattern
file_list2 = [file for file in file_list2 if not exclude_pattern.search(file)]

# Check if any files are left after exclusion
if not file_list2:
    raise ValueError("No files found after excluding the specified pattern.")


# Load dimensions from the first file
first_file2 = file_list2[0]
with Dataset(first_file2, 'r') as nc2:
    # Get dimensions and shape from pm25_tot
    pm25_tot_dims2 = nc.variables['pm25_ave'].dimensions
    dimensions2 = {dim: len(nc2.dimensions[dim]) for dim in pm25_tot_dims2}
    dim_names2 = list(dimensions2.keys())

    # Print each dimension's name and size
    print("Dimensions and sizes:")
    for dim_name2 in dim_names2:
        dim_size2 = dimensions[dim_name2]
        print(f"{dim_name2}: {dim_size2}")


# Initialize lists to collect data
data_pm25_tot = []
data_pm25_ave = []

# Load data from each file
for file in file_list:
    with Dataset(file, 'r') as nc:
        #if 'pm25_ave' in nc.variables and 'pm25_ave' in nc2.variables:
        if 'pm25_ave' in nc.variables:
            data_pm25_tot.append(nc.variables['pm25_ave'][:])
            #data_pm25_ave.append(nc2.variables['pm25_ave'][:])
            print(f"Find pm25_ave in file {file}")
        else:
            print(f"Missing pm25_ave in file {file}")

for file in file_list2:
    with Dataset(file, 'r') as nc2:
        if 'pm25_ave' in nc2.variables:
            #data_pm25_tot.append(nc.variables['pm25_ave'][:])
            data_pm25_ave.append(nc2.variables['pm25_ave'][:])
            print(f"Find pm25_ave in file {file}")
        else:
            print(f"Missing pm25_ave in file {file}")


# Convert lists to numpy arrays for easier manipulation
data_pm25_tot = np.concatenate(data_pm25_tot, axis=0)  # Shape: (n_files, ... ) depending on dimensions
data_pm25_ave = np.concatenate(data_pm25_ave, axis=0)  # Shape: (n_files, ... ) depending on dimensions

# Calculate statistics
mean_pm25_tot = np.mean(data_pm25_tot, axis=0)  # Shape: (dimensions excluding time)
std_pm25_tot = np.std(data_pm25_tot, axis=0)
mean_pm25_ave = np.mean(data_pm25_ave, axis=0)
std_pm25_ave = np.std(data_pm25_ave, axis=0)

# Calculate the mean and std of the difference between pm25_tot and pm25_ave
difference = data_pm25_ave - data_pm25_tot
mean_diff = np.mean(difference, axis=0)
std_diff = np.std(difference, axis=0)

# Write results to a new NetCDF file
output_file = 'output_statistics.nc'
with Dataset(output_file, 'w', format='NETCDF4') as nc_out:
    # Create dimensions
    for dim_name in dim_names:
        nc_out.createDimension(dim_name, dimensions[dim_name])

    # Create variables
    nc_out.createVariable('mean_pm25_da', np.float32, dim_names)
    nc_out.createVariable('std_pm25_da', np.float32, dim_names)
    nc_out.createVariable('mean_pm25_con', np.float32, dim_names)
    nc_out.createVariable('std_pm25_con', np.float32, dim_names)
    nc_out.createVariable('mean_diff', np.float32, dim_names)
    nc_out.createVariable('std_diff', np.float32, dim_names)

    # Write data to variables
    nc_out.variables['mean_pm25_da'][:] = mean_pm25_tot
    nc_out.variables['std_pm25_da'][:] = std_pm25_tot
    nc_out.variables['mean_pm25_con'][:] = mean_pm25_ave
    nc_out.variables['std_pm25_con'][:] = std_pm25_ave
    nc_out.variables['mean_diff'][:] = mean_diff
    nc_out.variables['std_diff'][:] = std_diff

print("Statistics written to output_statistics.nc")


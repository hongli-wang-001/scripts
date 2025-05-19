import numpy as np
import glob
from netCDF4 import Dataset

# Find all matching files
file_list = sorted(glob.glob('aod.2020*.nc'))

# Loop through each file and compute stats
for file in file_list:
    with Dataset(file, 'r') as ds:
        aod = ds.variables['aod550nm'][:]
        aod_flat = aod.flatten()
        aod_valid = aod_flat[~np.isnan(aod_flat)]  # Exclude NaNs
        
        mean_val = np.mean(aod_valid)
        std_val = np.std(aod_valid)

        print(f"{file}: mean = {mean_val:.4f}, std = {std_val:.4f}, count = {len(aod_valid)}")


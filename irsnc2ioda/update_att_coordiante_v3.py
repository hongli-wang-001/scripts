from netCDF4 import Dataset
import numpy as np
import shutil
import os

old_fname = "derived_rad_v15.ioda.nc"
new_fname = "derived_rad_v15_fixed_channel.ioda.nc"

# Copy file first
if os.path.exists(new_fname):
    os.remove(new_fname)
shutil.copy(old_fname, new_fname)

# Open copy in read/write
nc = Dataset(new_fname, mode='r+')

# Find Channel dimension in root
if 'Channel' not in nc.dimensions:
    raise KeyError("Dimension 'Channel' not found in root")

nch = nc.dimensions['Channel'].size
print("Channel dimension size:", nch)

# Check variable
if 'Channel' not in nc.variables:
    # Create it
    ch_var = nc.createVariable('Channel', np.int32, ('Channel',))
else:
    ch_var = nc.variables['Channel']

# Fill the values 1..nch
channel_numbers = np.arange(1, nch+1, dtype=np.int32)
ch_var[:] = channel_numbers

# Optional attributes
ch_var.long_name = "Channel index"
ch_var.units = ""
if '_FillValue' in ch_var.ncattrs():
    ch_var._FillValue = ch_var.getncattr('_FillValue')

nc.close()

print("Wrote new file:", new_fname)


import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np

# Define the NetCDF file and variable names
netcdf_file = 'file1.nc'
var1_name = 'aecj'
var2_name = 'aorgcj'

# Open the NetCDF file
with nc.Dataset(netcdf_file, 'r') as ds:
    # Read the variables
    var1 = ds.variables[var1_name][:]
    var2 = ds.variables[var2_name][:]
    
    # Check if dimensions match
    if var1.shape != var2.shape:
        raise ValueError("The dimensions of the variables do not match.")

    # Flatten the variables if they are multidimensional
    var1_flat = var1.flatten()
    var2_flat = var2.flatten()

    # Remove NaNs and infinities
    mask = np.isfinite(var1_flat) & np.isfinite(var2_flat)
    var1_clean = var1_flat[mask]
    var2_clean = var2_flat[mask]

# Create the scatter plot
plt.figure(figsize=(8, 6))
plt.scatter(var1_clean, var2_clean, alpha=0.5, edgecolors='w', s=20)
plt.xlabel(var1_name)
plt.ylabel(var2_name)
plt.title(f'Scatter Plot of {var1_name} vs. {var2_name}')
plt.grid(True)

# Save the plot to a file
plt.savefig('scatter_plot.png')  # Change the file name and extension as needed

print("Plot saved as 'scatter_plot.png'")


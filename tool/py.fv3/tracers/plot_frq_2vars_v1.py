import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Define the NetCDF file and variable names
netcdf_file = 'file1.nc'
variable_names = ['ano3k', 'anh4k']

# Open the NetCDF file and read the variables
with nc.Dataset(netcdf_file, 'r') as ds:
    # Read the 4D arrays from the specified variables
    data_ano3k = ds.variables[variable_names[0]][:]
    data_anh4k = ds.variables[variable_names[1]][:]
    
    # Ensure the data are 4D
    if data_ano3k.ndim != 4 or data_anh4k.ndim != 4:
        raise ValueError("One or both variables do not contain a 4D array.")
    
    # Check if the time dimension size is 1
    if data_ano3k.shape[0] != 1 or data_anh4k.shape[0] != 1:
        raise ValueError("The time dimension size is not 1.")
    
    # Extract the 3D arrays for time dimension 1
    data_ano3k_3d = data_ano3k[0, :, :, :]  # Shape: (zaxis_1, yaxis_1, xaxis_1)
    data_anh4k_3d = data_anh4k[0, :, :, :]  # Shape: (zaxis_1, yaxis_1, xaxis_1)

    # Flatten the data arrays to 1D
    data_ano3k_flat = data_ano3k_3d.flatten()
    data_anh4k_flat = data_anh4k_3d.flatten()

# Plot the histograms
plt.figure(figsize=(12, 8))

# Histogram for ano3k
plt.hist(data_ano3k_flat, bins=50, alpha=0.5, label='ano3k', color='blue')

# Histogram for anh4k
plt.hist(data_anh4k_flat, bins=50, alpha=0.5, label='anh4k', color='red')

plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Frequency Distribution of ano3k and anh4k')
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig('frequency_distribution_plot.png')  # Save the figure with the desired filename
print("Figure saved as 'frequency_distribution_plot.png'.")

# Optionally show the plot (can be omitted if running in a non-interactive environment)
plt.show()


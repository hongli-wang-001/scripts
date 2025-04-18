import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Define the NetCDF file and variable name
netcdf_file = 'file1.nc'
variable_name = 'asoil'

# Open the NetCDF file and read the variable
with nc.Dataset(netcdf_file, 'r') as ds:
    # Read the 4D array from the specified variable
    data = ds.variables[variable_name][:]
    
    # Ensure the data is 4D
    if data.ndim != 4:
        raise ValueError("The variable does not contain a 4D array.")
    
    # Check if the time dimension size is 1
    if data.shape[0] != 1:
        raise ValueError("The time dimension size is not 1.")
    
    # Extract the 3D array for time dimension 1
    data_3d = data[0, :, :, :]  # Shape: (zaxis_1, yaxis_1, xaxis_1)

    # Calculate median over xaxis_1 and yaxis_1 dimensions
    median_data = np.median(data_3d, axis=(1, 2))  # Shape: (zaxis_1)

# Plot the median data
z_dim = median_data.shape[0]
z = np.arange(z_dim)

plt.figure(figsize=(10, 6))
plt.plot(z, median_data, marker='o', linestyle='-', color='r')
plt.xlabel('Z Dimension')
plt.ylabel('Median Value')
plt.title('Median Data Over xaxis_1 and yaxis_1 Dimensions')
plt.grid(True)

# Save the plot to a file
plt.savefig('median_data_plot.png')  # Save the figure with the desired filename
print("Figure saved as 'median_data_plot.png'.")

# Optionally show the plot (can be omitted if running in a non-interactive environment)
plt.show()


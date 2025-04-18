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

    # Calculate mean and median over xaxis_1 and yaxis_1 dimensions
    mean_data = np.mean(data_3d, axis=(1, 2))  # Shape: (zaxis_1)
    median_data = np.median(data_3d, axis=(1, 2))  # Shape: (zaxis_1)

# Plot the data with z dimension on vertical axis
z_dim = mean_data.shape[0]
z = np.arange(z_dim)

plt.figure(figsize=(12, 8))
plt.plot(mean_data, z, marker='o', linestyle='-', color='b', label='Average')
plt.plot(median_data, z, marker='x', linestyle='--', color='r', label='Median')
plt.ylabel('Z Dimension')
plt.xlabel('Value')
plt.title('Average and Median Data Over xaxis_1 and yaxis_1 Dimensions')
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig('average_median_data_plot_vertical_z.png')  # Save the figure with the desired filename
print("Figure saved as 'average_median_data_plot_vertical_z.png'.")

# Optionally show the plot (can be omitted if running in a non-interactive environment)
plt.show()


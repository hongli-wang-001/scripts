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

    # Calculate mean and standard deviation over xaxis_1 and yaxis_1 dimensions
    mean_ano3k = np.mean(data_ano3k_3d, axis=(1, 2))  # Shape: (zaxis_1)
    std_ano3k = np.std(data_ano3k_3d, axis=(1, 2))  # Shape: (zaxis_1)
    
    mean_anh4k = np.mean(data_anh4k_3d, axis=(1, 2))  # Shape: (zaxis_1)
    std_anh4k = np.std(data_anh4k_3d, axis=(1, 2))  # Shape: (zaxis_1)

# Reverse the Z dimension
z_dim = mean_ano3k.shape[0]
z = np.arange(z_dim)  # Original z dimension order

# Plot the data
plt.figure(figsize=(12, 8))

# Plot average and standard deviation for ano3k
plt.errorbar(mean_ano3k, z, xerr=std_ano3k, label='ano3k', fmt='-o', color='blue', capsize=5)

# Plot average and standard deviation for anh4k
plt.errorbar(mean_anh4k, z, xerr=std_anh4k, label='anh4k', fmt='-x', color='red', capsize=5)

# Add labels and title
plt.xlabel('Value')
plt.ylabel('Z Dimension')
plt.title('Average and Standard Deviation of ano3k and anh4k Over xaxis_1 and yaxis_1 Dimensions')
plt.legend()
plt.grid(True)

# Invert the y-axis to have larger Z values at the bottom
plt.gca().invert_yaxis()

# Annotate the plot with average and std values
for i in range(len(z)):
    plt.text(mean_ano3k[i], z[i], f'{mean_ano3k[i]:.2f}\n± {std_ano3k[i]:.2f}', 
             color='blue', ha='center', va='center', fontsize=8, bbox=dict(facecolor='white', edgecolor='blue', boxstyle='round,pad=0.3'))
    plt.text(mean_anh4k[i], z[i], f'{mean_anh4k[i]:.2f}\n± {std_anh4k[i]:.2f}', 
             color='red', ha='center', va='center', fontsize=8, bbox=dict(facecolor='white', edgecolor='red', boxstyle='round,pad=0.3'))

# Save the plot to a file
plt.savefig('average_with_error_bars_and_annotations_plot.png')  # Save the figure with the desired filename
print("Figure saved as 'average_with_error_bars_and_annotations_plot.png'.")

# Optionally show the plot (can be omitted if running in a non-interactive environment)
plt.show()


import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Define the NetCDF files and variable names
netcdf_files = ['file1.nc', 'file2.nc', 'file3.nc']
variable_names = ['var1', 'var2', 'var3']  # Replace with actual variable names

# Read the one-dimensional arrays from the NetCDF files
data = []
for file, var in zip(netcdf_files, variable_names):
    with nc.Dataset(file, 'r') as ds:
        array = ds.variables[var][:]
        # Ensure the array is one-dimensional
        if array.ndim != 1:
            raise ValueError(f"The variable {var} in file {file} is not one-dimensional.")
        data.append(array)

# Combine all data to determine a common bin range
all_data = np.concatenate(data)

# Determine the bins based on the combined data
num_bins = 30  # Number of bins
bins = np.linspace(np.min(all_data), np.max(all_data), num_bins + 1)

# Plot the histograms
plt.figure(figsize=(12, 8))

# Plot histogram for each variable
colors = ['blue', 'green', 'red']
for i, (array, color) in enumerate(zip(data, colors)):
    plt.hist(array, bins=bins, alpha=0.6, color=color, label=f'Variable {i+1}')

# Add labels and title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Frequency Distribution of Three Variables')
plt.legend()
plt.grid(True)

# Save the plot to a file
plt.savefig('frequency_distribution_plot.png')  # Save the figure with the desired filename
print("Figure saved as 'frequency_distribution_plot.png'.")

# Optionally show the plot (can be omitted if running in a non-interactive environment)
plt.show()


import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob
import os

# Define file and variable details
input_pattern = '../3dvar*.nc'  # Adjust this pattern to match your file path
output_file = 'frequency_info.nc'
variable_name = 'particulatematter2p5Surface'
groups = ['oman', 'ombg']
missing_value = -3.368795e+38  # Missing value placeholder
max_limit = 3000  # Maximum value limit

# Function to extract and mask data from a NetCDF file
def extract_data(nc_file, group_name, var_name, missing_value):
    with nc.Dataset(nc_file, 'r') as dataset:
        data = dataset.groups[group_name][var_name][:]
        data = np.ma.masked_equal(data, missing_value)  # Mask missing values
    return data

# Function to plot scatter plot and print statistics
def plot_scatter(x_data, y_data, title, filename):
    x_data_flat = x_data.compressed()  # Get array without masked values
    y_data_flat = y_data.compressed()  # Get array without masked values

    # Calculate statistics
    x_min = np.min(x_data_flat)
    x_max = np.min(max(np.max(x_data_flat), max_limit))  # Check max limit
    y_min = np.min(y_data_flat)
    y_max = np.min(max(np.max(y_data_flat), max_limit))  # Check max limit
    x_mean = np.mean(x_data_flat)
    x_std_dev = np.std(x_data_flat)
    y_mean = np.mean(y_data_flat)
    y_std_dev = np.std(y_data_flat)

    # Plot scatter
    plt.figure(figsize=(10, 6))
    plt.scatter(x_data_flat, y_data_flat, alpha=0.5, color='blue', label='Data')
    plt.xlabel('Oman PM2.5')
    plt.ylabel('Ombg PM2.5')
    plt.title(title)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.grid(True)

    # Add text for statistics
    plt.text(0.95, 0.95, f'X Min: {x_min:.2f}\nX Max: {x_max:.2f}\nX Mean: {x_mean:.2f}\nX Std Dev: {x_std_dev:.2f}\n\n'
                         f'Y Min: {y_min:.2f}\nY Max: {y_max:.2f}\nY Mean: {y_mean:.2f}\nY Std Dev: {y_std_dev:.2f}',
             horizontalalignment='right', verticalalignment='top',
             transform=plt.gca().transAxes,
             fontsize=12, color='black', bbox=dict(facecolor='white', alpha=0.5))

    plt.legend()

    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()

    # Print statistics to the console
    print(f"{title} Statistics:")
    print(f"X Min: {x_min:.2f}")
    print(f"X Max: {x_max:.2f}")
    print(f"X Mean: {x_mean:.2f}")
    print(f"X Std Dev: {x_std_dev:.2f}")
    print(f"Y Min: {y_min:.2f}")
    print(f"Y Max: {y_max:.2f}")
    print(f"Y Mean: {y_mean:.2f}")
    print(f"Y Std Dev: {y_std_dev:.2f}")

# Get list of input files
input_files = glob.glob(input_pattern)

# Prepare to aggregate data
data_dict = {group: [] for group in groups}

# Process each file
for file in input_files:
    for group in groups:
        data = extract_data(file, group, variable_name, missing_value)
        data_dict[group].append(data)

# Combine all data for each group
combined_data_dict = {group: np.concatenate(data_dict[group]) for group in groups}

# Plot scatter plots
plot_scatter(combined_data_dict['oman'], combined_data_dict['ombg'],
             'Scatter Plot: Oman vs Ombg PM2.5', 'oman_vs_ombg_scatter.png')

# Save statistics to new NetCDF file
with nc.Dataset(output_file, 'w', format='NETCDF4') as new_nc:
    # Create dimensions
    new_nc.createDimension('group', len(groups))

    # Create variables for statistics
    min_var = new_nc.createVariable('min_value', 'f4', ('group',))
    max_var = new_nc.createVariable('max_value', 'f4', ('group',))
 


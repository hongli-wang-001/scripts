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

# Function to extract data from a NetCDF file
def extract_data(nc_file, group_name, var_name):
    with nc.Dataset(nc_file, 'r') as dataset:
        data = dataset.groups[group_name][var_name][:]
    return data

# Function to plot scatter plot and print statistics
def plot_scatter(data, title, filename, color, label):
    data_flat = data.ravel()

    # Calculate statistics
    data_min = np.min(data_flat)
    data_max = np.max(data_flat)
    mean = np.mean(data_flat)
    std_dev = np.std(data_flat)

    # Plot scatter
    plt.figure(figsize=(10, 6))
    plt.scatter(range(len(data_flat)), data_flat, color=color, alpha=0.5, label=label)
    plt.xlabel('Index')
    plt.ylabel('Value')
    plt.title(title)
    plt.grid(True)

    # Add text for statistics
    plt.text(0.95, 0.95, f'Min: {data_min:.2f}\nMax: {data_max:.2f}\nMean: {mean:.2f}\nStd Dev: {std_dev:.2f}',
             horizontalalignment='right', verticalalignment='top',
             transform=plt.gca().transAxes,
             fontsize=12, color=color, bbox=dict(facecolor='white', alpha=0.5))

    plt.legend()
    
    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()

    # Print statistics to the console
    print(f"{title} Statistics:")
    print(f"Min: {data_min:.2f}")
    print(f"Max: {data_max:.2f}")
    print(f"Mean: {mean:.2f}")
    print(f"Standard Deviation: {std_dev:.2f}")

# Get list of input files
input_files = glob.glob(input_pattern)

# Prepare to aggregate data
data_dict = {group: [] for group in groups}

# Process each file
for file in input_files:
    for group in groups:
        data = extract_data(file, group, variable_name)
        data_dict[group].append(data)

        # Plot scatter plots for each file
        plot_scatter(data, f'Scatter Plot - {group} ({os.path.basename(file)})', f'{group}_{variable_name}_{os.path.basename(file)}_scatter.png', 'blue' if group == 'oman' else 'orange', group)

# Combine all data for each group and plot combined scatter plots
for group in groups:
    combined_data = np.concatenate(data_dict[group])
    color = 'blue' if group == 'oman' else 'orange'
    label = group
    plot_scatter(combined_data, f'Combined Scatter Plot - {group}', f'{group}_{variable_name}_combined_scatter.png', color, label)

# Save statistics to new NetCDF file
with nc.Dataset(output_file, 'w', format='NETCDF4') as new_nc:
    # Create dimensions
    new_nc.createDimension('group', len(groups))

    # Create variables for statistics
    min_var = new_nc.createVariable('min_value', 'f4', ('group',))
    max_var = new_nc.createVariable('max_value', 'f4', ('group',))
    mean_var = new_nc.createVariable('mean_value', 'f4', ('group',))
    std_dev_var = new_nc.createVariable('std_dev', 'f4', ('group',))

    # Write statistics to NetCDF
    for i, group in enumerate(groups):
        combined_data = np.concatenate(data_dict[group])
        min_var[i] = np.min(combined_data)
        max_var[i] = np.max(combined_data)
        mean_var[i] = np.mean(combined_data)
        std_dev_var[i] = np.std(combined_data)

    print(f'Statistics for both groups saved to {output_file}')

print('Processing completed.')


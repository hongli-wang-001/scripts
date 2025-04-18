import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob
import os

# Define file and variable details
input_pattern = '../3dvar_nicas_pm25_2020090*nc'  # Change to your file path and pattern
output_file = 'frequency_info.nc'
variable_name = 'particulatematter2p5Surface'
groups = ['oman', 'ombg']

# Define bin range and size
bin_range = (-60, 60)
bin_size = 40
bins = np.linspace(bin_range[0], bin_range[1], bin_size)

# Function to extract data
def extract_data(nc_file, group_name, var_name):
    with nc.Dataset(nc_file, 'r') as dataset:
        data = dataset.groups[group_name][var_name][:]
    return data

# Function to plot scatter plot and print statistics
def plot_scatter(data, bins, title, filename, color, label, position=None, width=0.35):
    # Flatten the data for scatter plot
    data_flat = data.ravel()

    mean = np.mean(data_flat)
    std_dev = np.std(data_flat)

    # Create scatter plot
    plt.figure(figsize=(12, 7))
    plt.scatter(data_flat, np.zeros_like(data_flat), color=color, alpha=0.7, label=label, edgecolor='k')

    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title(title)
    plt.grid(True)

    # Add text for statistics
    plt.text(0.95, 0.95, f'Mean: {mean:.2f}\nStd Dev: {std_dev:.2f}',
             horizontalalignment='right', verticalalignment='top',
             transform=plt.gca().transAxes,
             fontsize=12, color=color, bbox=dict(facecolor='white', alpha=0.5))

    plt.legend()
    
    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()

# Get list of input files
input_files = glob.glob(input_pattern)

# Prepare to aggregate data
data_dict = {group: [] for group in groups}

# Process each file
for file in input_files:
    for group in groups:
        data = extract_data(file, group, variable_name)
        data_dict[group].append(data)

        # Plot individual scatter plots
        plot_scatter(data, bins, f'Scatter Plot - {group} ({os.path.basename(file)})', f'{group}_{variable_name}_{os.path.basename(file)}_scatter.png', 'blue' if group == 'oman' else 'orange', group)

# Combine all data for each group and plot combined scatter plots
plt.figure(figsize=(12, 7))
for group in groups:
    combined_data = np.concatenate(data_dict[group])
    color = 'blue' if group == 'oman' else 'orange'
    label = group
    plot_scatter(combined_data, bins, 'Combined Scatter Plot', None, color, label)

plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Combined Scatter Plot')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('combined_scatter_plot.png', dpi=300)
plt.close()

# Save bin and frequency info to new NetCDF file
with nc.Dataset(output_file, 'w', format='NETCDF4') as new_nc:
    # Create dimensions
    new_nc.createDimension('bins', len(bins) - 1)
    new_nc.createDimension('group', len(groups))

    # Create variables for bins and frequencies
    bin_edges_var = new_nc.createVariable('bin_edges', 'f4', ('bins',))
    frequencies_var = new_nc.createVariable('frequencies', 'i4', ('bins', 'group'))

    # Compute frequencies for each group
    all_frequencies = []
    for group in groups:
        combined_data = np.concatenate(data_dict[group])
        frequencies, _ = np.histogram(combined_data.ravel(), bins=bins)
        all_frequencies.append(frequencies)

    # Write data to NetCDF
    bin_edges_var[:] = bins[:-1]
    frequencies_var[:, 0] = all_frequencies[0]  # Frequencies for 'oman'
    frequencies_var[:, 1] = all_frequencies[1]  # Frequencies for 'ombg'

    print(f'Bin and frequency information for both groups saved to {output_file}')

print('Processing completed.')


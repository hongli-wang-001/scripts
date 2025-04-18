import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# Define file and variable details
input_file = './3dvar_nicas_pm25.nc'
output_file = 'frequency_info.nc'
variable_name = 'particulatematter2p5Surface'
groups = ['oman', 'ombg']

# Function to extract data
def extract_data(nc_file, group_name, var_name):
    with nc.Dataset(nc_file, 'r') as dataset:
        data = dataset.groups[group_name][var_name][:]
    return data

# Function to plot histogram and print statistics
def plot_histogram(data, bins, title, filename, color, label, position=None,text_y_position=None, width=0.35):
    frequencies, bin_edges = np.histogram(data.ravel(), bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    mean = np.mean(data)
    std_dev = np.std(data)

    if position is not None:
        plt.bar(bin_centers + position, frequencies, width=width, color=color, alpha=0.7, label=label, edgecolor='k')
    else:
        plt.hist(data.ravel(), bins=bins, edgecolor='k', alpha=0.7, color=color, label=label)

    if text_y_position is not None:
       y_pos = text_y_position 
    else:
       y_pos = 0.95


    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.grid(True)

    # Print statistics
    print(f"{title} Statistics:")
    print(f"Mean: {mean:.2f}")
    print(f"Standard Deviation: {std_dev:.2f}")

    # Add text for statistics
    plt.text(0.25, y_pos, f'Mean: {mean:.2f}\nStd Dev: {std_dev:.2f}',
             horizontalalignment='right', verticalalignment='top',
             transform=plt.gca().transAxes,
             fontsize=12, color=color, bbox=dict(facecolor='white', alpha=0.5))

    plt.legend()
    
    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()

# Extract data and plot histograms
data_dict = {}
for group in groups:
    data = extract_data(input_file, group, variable_name)
    data_dict[group] = data

    # Define bins
    bin_range = (-60, 60)
    bin_size = 40
    bins = np.linspace(bin_range[0], bin_range[1], bin_size)

    # Plot individual histograms
    plot_histogram(data, bins, f'Frequency Distribution - {group}', f'{group}_{variable_name}_histogram.png', 'blue' if group == 'oman' else 'orange', group)

# Plot combined histogram
plt.figure(figsize=(12, 7))
width = 0.5  # Width of bars

# Define position offsets for each group
positions = {
    'oman': -width / 2,
    'ombg': width / 2
}

text_y_positions = {
    'oman': 0.95,
    'ombg': 0.85 
}

for group in groups:
    data = data_dict[group]
    color = 'blue' if group == 'oman' else 'orange'
    label = group
    # Calculate histogram data
    frequencies, bin_edges = np.histogram(data.ravel(), bins=bins)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    plot_histogram(data, bins, 'Combined Frequency Distribution', None, color, label, position=positions[group], text_y_position=text_y_positions[group], width=width)

plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Combined Frequency Distribution')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('combined_histogram.png', dpi=300)
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
        data = data_dict[group]
        frequencies, _ = np.histogram(data.ravel(), bins=bins)
        all_frequencies.append(frequencies)

    # Write data to NetCDF
    bin_edges_var[:] = bins[:-1]
    frequencies_var[:, 0] = all_frequencies[0]  # Frequencies for 'oman'
    frequencies_var[:, 1] = all_frequencies[1]  # Frequencies for 'ombg'

    print(f'Bin and frequency information for both groups saved to {output_file}')

print('Processing completed.')


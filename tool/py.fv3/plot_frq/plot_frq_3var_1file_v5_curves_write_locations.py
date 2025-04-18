import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Define the NetCDF file and groups with variable name
netcdf_file = 'file1.nc'
groups = ['ObsValue', 'hofx0', 'hofx1']  # Replace with actual group names
variable_name = 'particulatematter2p5Surface'  # Actual variable name

# Read the one-dimensional arrays from different groups in the NetCDF file
data = []
statistics = {}
bin_locations = {}
for group in groups:
    with nc.Dataset(netcdf_file, 'r') as ds:
        # Access the variable within the specified group
        array = ds.groups[group].variables[variable_name][:]
        # Ensure the array is one-dimensional
        if array.ndim != 1:
            raise ValueError(f"The variable {variable_name} in group {group} is not one-dimensional.")
        data.append(array)
        
        # Calculate and store statistics for the current variable
        min_val = np.min(array)
        max_val = np.max(array)
        mean_val = np.mean(array)
        std_val = np.std(array)
        
        statistics[group] = {
            'min': min_val,
            'max': max_val,
            'mean': mean_val,
            'std': std_val
        }

        # Determine bin edges and count frequencies
        max_all_data = np.max(array)
        bins = np.linspace(0, max_all_data, num=50)  # Adjusted number of bins to 50
        counts, edges = np.histogram(array, bins=bins)
        bin_indices = [np.where((array >= edges[i]) & (array < edges[i+1]))[0] for i in range(len(edges)-1)]
        
        bin_locations[group] = bin_indices

# Print bin locations for each group
for group in groups:
    print(f"\nBin locations for group {group}:")
    for i, (start, end) in enumerate(zip(bins[:-1], bins[1:])):
        indices = bin_locations[group][i]
        print(f"  Bin {i+1} ({start:.2f} - {end:.2f}): {len(indices)} values")
        if len(indices) > 0:
            print(f"    Indices: {indices}")

# Plot the histograms
plt.figure(figsize=(12, 8))

# Plot histogram for each group
colors = ['blue', 'green', 'red']
for i, (array, color) in enumerate(zip(data, colors)):
    plt.hist(array, bins=bins, alpha=0.6, color=color, label=f'Group {groups[i]}')

# Add labels and title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Frequency Distribution of Variable Across Different Groups')
plt.legend()
plt.grid(True)

# Annotate the figure with statistics
textstr = '\n'.join([
    f"{group}:\n  Min: {stats['min']:.2f}\n  Max: {stats['max']:.2f}\n  Mean: {stats['mean']:.2f}\n  Std Dev: {stats['std']:.2f}"
    for group, stats in statistics.items()
])

# Add text annotation to the top-right corner
plt.gca().annotate(textstr, xy=(1.0, 1.0), xycoords='axes fraction',
                    fontsize=10, verticalalignment='top',
                    horizontalalignment='right', bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

# Save the histogram plot to a file
plt.savefig('frequency_distribution_histogram.png')  # Save the histogram figure with the desired filename
print("Histogram figure saved as 'frequency_distribution_histogram.png'.")

# Create frequency distribution curves
plt.figure(figsize=(12, 8))

# Plot frequency distribution for each group using curves
for i, (array, color) in enumerate(zip(data, colors)):
    # Compute histogram and density
    counts, edges = np.histogram(array, bins=bins)
    bin_centers = 0.5 * (edges[1:] + edges[:-1])
    # Plot curve
    plt.plot(bin_centers, counts, color=color, label=f'Group {groups[i]}')

# Add labels and title
plt.xlabel('Value')
plt.ylabel('Frequency')
plt.title('Frequency Distribution Curves of Variable Across Different Groups')
plt.legend()
plt.grid(True)

# Annotate the figure with statistics
textstr = '\n'.join([
    f"{group}:\n  Mean: {stats['mean']:.2f}\n  Std Dev: {stats['std']:.2f}"
    for group, stats in statistics.items()
])

# Add text annotation to the top-right corner
plt.gca().annotate(textstr, xy=(1.0, 1.0), xycoords='axes fraction',
                    fontsize=10, verticalalignment='top',
                    horizontalalignment='right', bbox=dict(boxstyle='round,pad=0.3', edgecolor='black', facecolor='white'))

# Save the KDE plot to a file
plt.savefig('frequency_distribution_curves.png')  # Save the curves figure with the desired filename
print("Curves figure saved as 'frequency_distribution_curves.png'.")

# Optionally show the plots (can be omitted if running in a non-interactive environment)
plt.show()


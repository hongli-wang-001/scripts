import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob

# Define file and variable details
input_file = './3dvar_nicas_pm25_to_2020091106.nc'  # Adjust this to your file path
output_file = 'frequency_info.nc'
variable_name = 'particulatematter2p5Surface'
groups = ['oman', 'ombg']
missing_value = -3.368795e+38  # Missing value placeholder
max_limit = 3000  # Maximum value limit
min_limit = -3000  # Minimum value limit

# Function to extract and mask data from a NetCDF file
def extract_data(nc_file, group_name, var_name, missing_value):
    with nc.Dataset(nc_file, 'r') as dataset:
        data = dataset.groups[group_name][var_name][:]
        data_masked = np.ma.masked_equal(data, missing_value)
        data_masked = np.ma.clip(data_masked, min_limit, max_limit)
    return data_masked

# Function to print statistics
def print_statistics(data, group_name):
    data_flat = data.compressed()
    filtered_data_flat = data_flat[np.abs(data_flat) > 100]
    print(f"Filtered data for {group_name}:")
    print(filtered_data_flat)
    if len(data_flat) == 0:
        print(f"Statistics for {group_name}: No valid data available.")
        return

    min_val = np.min(data_flat)
    max_val = np.max(data_flat)

    min_val = max(min_val, min_limit)
    max_val = min(max_val, max_limit)

    mean_val = np.mean(data_flat)
    std_dev_val = np.std(data_flat)

    mean_val = float('nan') if np.isinf(mean_val) else mean_val
    std_dev_val = float('nan') if np.isinf(std_dev_val) else std_dev_val

    print(f"Statistics for {group_name}:")
    print(f"Min: {min_val:.2f}")
    print(f"Max: {max_val:.2f}")
    print(f"Mean: {mean_val:.2f}" if not np.isnan(mean_val) else "Mean: NaN")
    print(f"Standard Deviation: {std_dev_val:.2f}" if not np.isnan(std_dev_val) else "Standard Deviation: NaN")
    print()

# Function to plot scatter plot with density
def plot_scatter(x_data, y_data, title, filename):
    x_data_flat = x_data.compressed()
    y_data_flat = y_data.compressed()

    if len(x_data_flat) == 0 or len(y_data_flat) == 0:
        print(f"Scatter plot cannot be generated. No valid data available for {title}.")
        return

    # Filtering data for absolute values greater than 100
    filtered_x_data_flat = x_data_flat[np.abs(x_data_flat) > 100]
    filtered_y_data_flat = y_data_flat[np.abs(y_data_flat) > 100]

    print(f"x_data_flat_ABS > 100")
    print(filtered_x_data_flat)
    print(f"y_data_flat_ABS > 100")
    print(filtered_y_data_flat)

    x_min = np.min(x_data_flat)
    x_max = np.max(x_data_flat)
    y_min = np.min(y_data_flat)
    y_max = np.max(y_data_flat)

    x_min = max(x_min, min_limit)
    x_max = min(x_max, max_limit)
    y_min = max(y_min, min_limit)
    y_max = min(y_max, max_limit)

    x_mean = np.mean(x_data_flat)
    x_std_dev = np.std(x_data_flat)
    y_mean = np.mean(y_data_flat)
    y_std_dev = np.std(y_data_flat)

    x_mean = float('nan') if np.isinf(x_mean) else x_mean
    x_std_dev = float('nan') if np.isinf(x_std_dev) else x_std_dev
    y_mean = float('nan') if np.isinf(y_mean) else y_mean
    y_std_dev = float('nan') if np.isinf(y_std_dev) else y_std_dev

    plt.figure(figsize=(10, 6))
    hexbin = plt.hexbin(x_data_flat, y_data_flat, gridsize=50, cmap='Blues', mincnt=1, alpha=0.75)
#    plt.colorbar(label='Density')
# Custom scaling (e.g., top 90% density)
    max_density = np.percentile(hexbin.get_array(), 90)
    plt.colorbar(label='Density', extend='max')
# Customize the color limits
    plt.clim(0, max_density)  # Set color limits based on the percentile

    plt.xlabel('Oman PM2.5')
    plt.ylabel('Ombg PM2.5')
    plt.title(title)
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.grid(True)

    # Add statistics text
    plt.text(0.95, 0.95, f'X Min: {x_min:.2f}\nX Max: {x_max:.2f}\nX Mean: {x_mean:.2f}\nX Std Dev: {x_std_dev:.2f}\n\n'
                         f'Y Min: {y_min:.2f}\nY Max: {y_max:.2f}\nY Mean: {y_mean:.2f}\nY Std Dev: {y_std_dev:.2f}',
             horizontalalignment='right', verticalalignment='top',
             transform=plt.gca().transAxes,
             fontsize=12, color='black', bbox=dict(facecolor='white', alpha=0.5))

    plt.legend(['Data'], loc='best')

    if filename:
        plt.savefig(filename, dpi=300)
        plt.close()

    print(f"{title} Statistics:")
    print(f"X Min: {x_min:.2f}")
    print(f"X Max: {x_max:.2f}")
    print(f"X Mean: {x_mean:.2f}")
    print(f"X Std Dev: {x_std_dev:.2f}")
    print(f"Y Min: {y_min:.2f}")
    print(f"Y Max: {y_max:.2f}")
    print(f"Y Mean: {y_mean:.2f}")
    print(f"Y Std Dev: {y_std_dev:.2f}")

# Process the single file
data_dict = {group: [] for group in groups}
missing_values_count = {group: 0 for group in groups}

for group in groups:
    data = extract_data(input_file, group, variable_name, missing_value)
    print(f"data_ABS > 100")
    filtered_data = data[np.abs(data) > 100]
    print(filtered_data)
    print_statistics(data, group)
    data_dict[group] = data
    missing_count = np.sum(data.mask)
    missing_values_count[group] += missing_count
    print(f"Number of missing values for {group} in file {input_file}: {missing_count}")

# Extract data for plotting
x_data = data_dict['oman']
y_data = data_dict['ombg']

# Generate scatter plot with density
plot_scatter(x_data, y_data, "Scatter Plot of Oman vs Ombg PM2.5", "scatter_plot_to_2020091106.png")

# Save statistics to new NetCDF file
with nc.Dataset(output_file, 'w', format='NETCDF4') as new_nc:
    # Create dimensions
    new_nc.createDimension('group', len(groups))

    # Create variables for statistics
    min_var = new_nc.createVariable('min_value', 'f4', ('group',))
    max_var = new_nc.createVariable('max_value', 'f4', ('group',))
    mean_var = new_nc.createVariable('mean_value', 'f4', ('group',))
    std_dev_var = new_nc.createVariable('std_dev_value', 'f4', ('group',))

    # Assign values
    min_var[:] = [np.min(data_dict[group].compressed()) for group in groups]
    max_var[:] = [np.max(data_dict[group].compressed()) for group in groups]
    mean_var[:] = [np.mean(data_dict[group].compressed()) for group in groups]
    std_dev_var[:] = [np.std(data_dict[group].compressed()) for group in groups]

    # Add attributes
    min_var.units = 'units'
    max_var.units = 'units'
    mean_var.units = 'units'
    std_dev_var.units = 'units'
    min_var.description = 'Minimum value of the data'
    max_var.description = 'Maximum value of the data'
    mean_var.description = 'Mean value of the data'
    std_dev_var.description = 'Standard deviation of the data'

print("Statistics saved to", output_file)


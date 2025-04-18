import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import glob

# Define file and variable details
input_pattern = '../3dvar_nicas_pm25_2020090412.nc'  # Adjust this pattern to match your file path
output_file = 'frequency_info.nc'
variable_name = 'particulatematter2p5Surface'
groups = ['oman', 'ombg']
missing_value = -3.368795e+38  # Missing value placeholder
max_limit = 3000  # Maximum value limit
default_min_limit = -3000  # Default minimum value limit

# Function to extract and mask data from a NetCDF file
def extract_data(nc_file, group_name, var_name, missing_value):
    with nc.Dataset(nc_file, 'r') as dataset:
        data = dataset.groups[group_name][var_name][:]
        print(f"Raw data from {group_name}:")
        # Print all data where absolute value > 100
        filtered_raw_data = data[np.abs(data) > 100]
        print(filtered_raw_data)
        
        data_masked = np.ma.masked_equal(data, missing_value)  # Mask missing values
        print(f"Masked data from {group_name}:")
        # Print all masked data where absolute value > 100
        filtered_masked_data = data_masked[np.abs(data_masked) > 100]
        print(filtered_masked_data)
        
    return data_masked

# Function to print statistics
def print_statistics(data, group_name):
    data_flat = data.compressed()  # Get array without masked values
    print(f"Compressed data for {group_name}:")
    # Print all data where absolute value > 100
    filtered_data_flat = data_flat[np.abs(data_flat) > 100]
    print(filtered_data_flat)

    if len(filtered_data_flat) == 0:  # Check if there are no valid data points
        print(f"Statistics for {group_name}: No valid data available.")
        return

    min_val = np.min(filtered_data_flat)
    max_val = np.max(filtered_data_flat)
    
    # Use min and max values from data for clipping
    min_val = max(min_val, default_min_limit)
    max_val = min(max_val, max_limit)
    
    mean_val = np.mean(filtered_data_flat)
    std_dev_val = np.std(filtered_data_flat)
    
    # Check for infinities and NaNs
    mean_val = float('nan') if np.isinf(mean_val) else mean_val
    std_dev_val = float('nan') if np.isinf(std_dev_val) else std_dev_val

    # Print valid statistics
    print(f"Statistics for {group_name}:")
    print(f"Min: {min_val:.2f}")
    print(f"Max: {max_val:.2f}")
    print(f"Mean: {mean_val:.2f}" if not np.isnan(mean_val) else "Mean: NaN")
    print(f"Standard Deviation: {std_dev_val:.2f}" if not np.isnan(std_dev_val) else "Standard Deviation: NaN")
    print()

# Function to plot scatter plot and print statistics
def plot_scatter(x_data, y_data, title, filename):
    x_data_flat = x_data.compressed()  # Get array without masked values
    y_data_flat = y_data.compressed()  # Get array without masked values

    if len(x_data_flat) == 0 or len(y_data_flat) == 0:  # Check if there's no valid data for plotting
        print(f"Scatter plot cannot be generated. No valid data available for {title}.")
        return

    # Calculate statistics
    x_min = np.min(x_data_flat)
    x_max = np.max(x_data_flat)
    y_min = np.min(y_data_flat)
    y_max = np.max(y_data_flat)
    
    # Use min and max values from data for clipping
    x_min = max(x_min, default_min_limit)
    x_max = min(x_max, max_limit)
    y_min = max(y_min, default_min_limit)
    y_max = min(y_max, max_limit)
    
    x_mean = np.mean(x_data_flat)
    x_std_dev = np.std(x_data_flat)
    y_mean = np.mean(y_data_flat)
    y_std_dev = np.std(y_data_flat)
    
    # Check for infinities and NaNs
    x_mean = float('nan') if np.isinf(x_mean) else x_mean
    x_std_dev = float('nan') if np.isinf(x_std_dev) else x_std_dev
    y_mean = float('nan') if np.isinf(y_mean) else y_mean
    y_std_dev = float('nan') if np.isinf(y_std_dev) else y_std_dev

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
missing_values_count = {group: 0 for group in groups}

# Process each file
for file in input_files:
    for group in groups:
        data = extract_data(file, group, variable_name, missing_value)
        data_dict[group].append(data)
        # Print statistics for each group
        print_statistics(data, group)
        # Count and print missing values
        missing_count = np.sum(data.mask)
        missing_values_count[group] += missing_count
        print(f"Number of missing values for {group} in file {file}: {missing_count}")

# Combine all data for each group
combined_data_dict = {group: np.concatenate(data_dict[group]) for group in groups}

# Print the number of missing values
for group, count in missing_values_count.items():
    print(f"Total Missing Values for {group}: {count}")

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
    mean_var = new_nc.createVariable('mean_value', 'f4', ('group',))
    std_dev_var = new_nc.createVariable('std_dev', 'f4', ('group',))
    missing_var = new_nc.createVariable('missing_values', 'i4', ('group',))

    # Write statistics to NetCDF
    for i, group in enumerate(groups):
        combined_data = combined_data_dict[group]
        min_val = np.min(combined_data)
        max_val = np.max(combined_data)
        min_val = max(min_val, default_min_limit)  # Apply min limit
        max_val = min(max_val, max_limit)  # Apply max limit
        mean_val = np.mean(combined_data)
        std_dev_val = np.std(combined_data)
        missing_count = missing_values_count[group]
        
        min_var[i] = min_val
        max_var[i] = max_val
        mean_var[i] = mean_val
        std_dev_var[i] = std_dev_val
        missing_var[i] = missing_count


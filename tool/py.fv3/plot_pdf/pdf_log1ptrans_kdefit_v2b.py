import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.stats import gaussian_kde

# Define file and variable details
netcdf_file = 'file1.nc'
groups = ['ObsValue','hofx0','hofx1']  # List of groups in the NetCDF file
variable_name = 'particulatematter2p5Surface'  # Variable name to extract

# Function to extract data from NetCDF
def extract_data(nc_file, group_name, var_name):
    """Extract data from a specified NetCDF file, group, and variable."""
    with Dataset(nc_file, 'r') as nc:
        data = nc.groups[group_name][var_name][:]
    return data

# Extract data for each group
data_dict = {}
for group in groups:
    data = extract_data(netcdf_file, group, variable_name)
    data_dict[group] = data

# Apply log transformation, handling zeros by adding a small constant
transformed_data_dict = {}
for group, data in data_dict.items():
    # Replace zeros with NaN to avoid log(0) issues
    data_cleaned = np.where(data > 0, data, np.nan)
    # Apply log1p for a more numerically stable transformation (log(1 + x))
    transformed_data_dict[group] = np.log1p(data_cleaned)

# Determine the bin range for log-transformed data
all_log_values = np.concatenate([transformed_data_dict[group].flatten() for group in groups])
all_log_values = all_log_values[~np.isnan(all_log_values)]  # Remove NaNs for fitting

if len(all_log_values) == 0:
    raise ValueError("No valid data points available for histogram and KDE fitting.")

# Fit Kernel Density Estimate
kde = gaussian_kde(all_log_values, bw_method='scott')  # 'scott' or 'silverman' for bandwidth selection

# Generate data for plotting the KDE fit
x = np.linspace(np.min(all_log_values), np.max(all_log_values), 1000)
pdf = kde(x)

# Determine bin range for plotting
bin_range_log = (np.min(all_log_values), np.max(all_log_values))
num_bins_log = 40
bins_log = np.linspace(bin_range_log[0], bin_range_log[1], num_bins_log + 1)

# Compute original bin edges
bins_original = np.expm1(bins_log)  # Inverse of log1p

# Plot frequency distribution with KDE fit
fig, axes = plt.subplots(nrows=1, ncols=len(groups), figsize=(15, 5), sharex='col', sharey='col')

# If only one subplot, axes will be a single Axes object
if len(groups) == 1:
    axes = [axes]

for i, group in enumerate(groups):
    ax = axes[i]
    # Plot histogram of log-transformed data
    counts, bins_log = np.histogram(transformed_data_dict[group].flatten(), bins=bins_log)
    bin_centers_log = (bins_log[:-1] + bins_log[1:]) / 2
    ax.bar(bin_centers_log, counts, width=np.diff(bins_log), alpha=0.7, edgecolor='black', label='Histogram')

    # Plot KDE fit
    ax.plot(x, pdf * np.max(counts) / np.max(pdf), 'r-', lw=2, label='KDE Fit')

    # Set axis labels and title
    ax.set_xlabel('Log(Value)')
    ax.set_ylabel('Frequency')
    #ax.set_title(f'{group} (Log Transformed)')
    ax.set_title(f'{group}')
    ax.legend()

    # Specify the maximum value for y-axis
    # Commented the two lines if you want automotic maximu of y-axis
    max_y = 150  # Set this to your desired maximum y-axis value
    ax.set_ylim(0, max_y)

    # Adding secondary x-axis for original values
    ax2 = ax.twiny()
    ax2.set_xlim(np.log1p(bin_range_log[0]), np.log1p(bin_range_log[1]))

    # Calculate tick positions and labels
    ticks_positions = np.log1p(bins_original)
    labels = [f'{int(val):,}' for val in bins_original]

    # Set ticks and labels, only every third label
    ax2.set_xticks(ticks_positions)
    ax2.set_xticklabels(labels)

    # Only show every third label
    num_labels = len(labels)
    ax2.set_xticks([ticks_positions[i] for i in range(0, num_labels, 6)])
    ax2.set_xticklabels([labels[i] for i in range(0, num_labels, 6)])

    ax2.set_xlabel('Original Value')

plt.tight_layout()
plt.savefig('frequency_distribution_kde_fit.png')
plt.show()


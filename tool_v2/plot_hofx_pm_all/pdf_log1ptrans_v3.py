import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.stats import norm

# Define file and variable details
netcdf_file = 'file1.nc'
groups = ['ObsValue']  # List of groups in the NetCDF file
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
    raise ValueError("No valid data points available for histogram and GMM fitting.")

# Fit a Gaussian Mixture Model manually using scipy
num_components = 2  # Adjust the number of components as needed
means = []
stds = []
weights = []

# Fit Gaussian distributions manually
for i in range(num_components):
    mean = np.mean(all_log_values)  # Simple mean as initial guess
    std = np.std(all_log_values)    # Simple std dev as initial guess
    weights.append(1 / num_components)
    means.append(mean)
    stds.append(std)

# Generate data for plotting the Gaussian fits
x = np.linspace(np.min(all_log_values), np.max(all_log_values), 1000)
pdf_gaussian = np.zeros_like(x)

for mean, std, weight in zip(means, stds, weights):
    pdf_gaussian += weight * norm.pdf(x, loc=mean, scale=std)

# Determine bin range for plotting
bin_range_log = (np.min(all_log_values), np.max(all_log_values))
num_bins_log = 40
bins_log = np.linspace(bin_range_log[0], bin_range_log[1], num_bins_log + 1)

# Plot frequency distribution with Gaussian fits
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

    # Plot Gaussian fits
    ax.plot(x, pdf_gaussian * np.max(counts) / np.max(pdf_gaussian), 'r-', lw=2, label='Gaussian Fit')

    # Set axis labels and title
    ax.set_xlabel('Log(Value)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{group} (Log Transformed)')
    ax.legend()

    # Annotate Gaussian parameters
    for j in range(num_components):
        ax.text(0.05, 0.9 - j * 0.1, f'Component {j + 1}:\nMean: {means[j]:.2f}\n'
                f'Std: {stds[j]:.2f}\nWeight: {weights[j]:.2f}',
                transform=ax.transAxes, fontsize=12, verticalalignment='top')

plt.tight_layout()
plt.savefig('frequency_distribution_gaussian_fit.png')
plt.show()


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

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

# Fit Gaussian Mixture Model
n_components = 3  # Number of Gaussian components (adjust as needed)
gmm = GaussianMixture(n_components=n_components, covariance_type='full')
gmm.fit(all_log_values.reshape(-1, 1))

# Get the parameters of the GMM
means = gmm.means_.flatten()
covariances = gmm.covariances_.flatten()
weights = gmm.weights_

# Generate data for plotting the Gaussian Mixture fit
x = np.linspace(np.min(all_log_values), np.max(all_log_values), 1000)
pdf = np.zeros_like(x)
for i in range(n_components):
    pdf += weights[i] * (1 / np.sqrt(2 * np.pi * covariances[i])) * np.exp(-0.5 * ((x - means[i]) ** 2) / covariances[i])

# Plot frequency distribution with Gaussian Mixture fit
fig, axes = plt.subplots(nrows=1, ncols=len(groups), figsize=(15, 5), sharex='col', sharey='col')

# If only one subplot, axes will be a single Axes object
if len(groups) == 1:
    axes = [axes]

for i, group in enumerate(groups):
    ax = axes[i]
    # Plot histogram of log-transformed data
    counts, bins_log = np.histogram(transformed_data_dict[group].flatten(), bins=40)
    ax.bar(bins_log[:-1], counts, width=np.diff(bins_log), alpha=0.7, edgecolor='black', label='Histogram')

    # Plot Gaussian Mixture fit
    bin_centers_log = (bins_log[:-1] + bins_log[1:]) / 2
    ax.plot(x, pdf * np.max(counts) / np.max(pdf), 'r-', lw=2, label='Gaussian Mixture Fit')

    # Set axis labels and title
    ax.set_xlabel('Log(Value)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{group} (Log Transformed)')

    ax.legend()

plt.tight_layout()
plt.savefig('frequency_distribution_gmm_fit.png')
plt.show()


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from scipy.stats import gaussian_kde

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

# Function to compute the Modified Z-Score
def modified_z_score(data, threshold=3.5):
    """Calculate Modified Z-Scores and return indices of outliers."""
    med = np.nanmedian(data)
    mad = np.nanmedian(np.abs(data - med))
    if mad == 0:
        return np.array([])
    z_scores = 0.6745 * (data - med) / mad
    return np.abs(z_scores) > threshold

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

# Identify outliers using Modified Z-Score
outlier_indices = modified_z_score(all_log_values)
cleaned_log_values = all_log_values[~outlier_indices]

if len(cleaned_log_values) == 0:
    raise ValueError("No valid data points remaining after removing outliers.")

# Fit Kernel Density Estimate on the cleaned data
kde_cleaned = gaussian_kde(cleaned_log_values, bw_method='scott')  # 'scott' or 'silverman' for bandwidth selection

# Generate data for plotting the KDE fit
x = np.linspace(np.min(all_log_values), np.max(all_log_values), 1000)
pdf = kde_cleaned(x)

# Determine bin range for plotting
bin_range_log = (np.min(all_log_values), np.max(all_log_values))
num_bins_log = 40
bins_log = np.linspace(bin_range_log[0], bin_range_log[1], num_bins_log + 1)

# Plot frequency distribution with KDE fit
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 5), sharex='col', sharey='col') 

# Ensure axes is a list if only one group
if len(groups) == 1:
    axes = [axes]

for i, group in enumerate(groups):
    # Original Data Plot
    ax = axes[i]
    ax.hist(transformed_data_dict[group].flatten(), bins=bins_log, alpha=0.7, edgecolor='black', label='Original Data')
    # KDE Fit on original data
    kde_original = gaussian_kde(all_log_values, bw_method='scott')
    pdf_original = kde_original(x)
    ax.plot(x, pdf_original * np.max(np.histogram(transformed_data_dict[group].flatten(), bins=bins_log)[0]) / np.max(pdf_original), 'r-', lw=2, label='KDE Fit (Original)')
    ax.set_xlabel('Log(Value)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{group} (Log Transformed)')
    ax.legend()

    # Cleaned Data Plot (without outliers)
    ax = axes[i+1]
    ax.hist(cleaned_log_values, bins=bins_log, alpha=0.7, edgecolor='black', label='Cleaned Data')
    # KDE Fit on cleaned data
    ax.plot(x, pdf * np.max(np.histogram(cleaned_log_values, bins=bins_log)[0]) / np.max(pdf), 'r-', lw=2, label='KDE Fit (Cleaned)')
    ax.set_xlabel('Log(Value)')
    ax.set_ylabel('Frequency')
    ax.set_title(f'{group} (Without Outliers)')
    ax.legend()

plt.tight_layout()
plt.savefig('frequency_distribution_with_and_without_outliers.png')
plt.show()

# Print parameters of KDE fits
print("KDE Parameters for Original Data:")
print(f"  Bandwidth method: scott")
print("KDE Parameters for Cleaned Data:")
p


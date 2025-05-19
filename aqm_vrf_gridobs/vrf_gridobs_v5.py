import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Open fcst.nc (read-only) and pm25.nc (append mode)
with nc.Dataset('fcst.nc', 'r') as fcst_ds, nc.Dataset('pm25.nc', 'a') as pm25_ds:
    # Read the specific slice from fcst.nc
    pm25_ave = fcst_ds.variables['pm25_ave'][0, 63, :, :]  # Shape: [225, 393]

    # Read pm25 from pm25.nc
    pm25 = pm25_ds.variables['pm25'][:, :]  # Shape: [225, 393]

    # Flatten and mask invalid (NaN or fill value -1)
    pred = pm25_ave.flatten()
    obs = pm25.flatten()

    valid_mask = (~np.isnan(pred)) & (~np.isnan(obs)) & (obs != -1) & (pred != -1)
    pred_valid = pred[valid_mask]
    obs_valid = obs[valid_mask]

    # Compute difference
    diff = pred_valid - obs_valid

    # Remove or mask invalid values if needed (e.g., NaNs)
    valid_diff = diff[~np.isnan(diff)]

    # Count valid data points
    num_valid_points = valid_mask.sum()
 
    # Calculate statistics
    mean_diff = np.mean(valid_diff)
    median_diff = np.median(valid_diff)
    std_diff = np.std(valid_diff)
    rms_diff = np.sqrt(np.mean(valid_diff ** 2))

    # Correlation (Pearson)
    if len(pred_valid) > 1:
        corr_matrix = np.corrcoef(pred_valid, obs_valid)
        correlation = corr_matrix[0, 1]
    else:
        correlation = np.nan

    # Index of Agreement (IOA)
    obs_mean = np.mean(obs_valid)
    numerator = np.sum((pred_valid - obs_valid) ** 2)
    denominator = np.sum((np.abs(pred_valid - obs_mean) + np.abs(obs_valid - obs_mean)) ** 2)
    ioa = 1 - (numerator / denominator) if denominator != 0 else np.nan

    # Print statistics
    print("Statistics for pm25_diff (pm25_ave[0,63,:,:] - pm25):")
    print(f"  Valid data points: {num_valid_points}")
    print(f"  Mean      : {mean_diff:.4f}")
    print(f"  Median    : {median_diff:.4f}")
    print(f"  Std Dev   : {std_diff:.4f}")
    print(f"  RMS       : {rms_diff:.4f}")
    print(f"  IOA       : {ioa:.4f}")
    print(f"  COR       : {correlation:.4f}")

    # Save the full difference field into pm25.nc
    full_diff = pm25_ave - pm25

    # If variable already exists, delete it to avoid conflict
    if 'pm25_diff' in pm25_ds.variables:
        diff_var = pm25_ds.variables['pm25_diff']
    else:
        diff_var = pm25_ds.createVariable('pm25_diff', 'f4', ('y', 'x'))

    diff_var[:, :] = full_diff 

    # Add metadata
    diff_var.units = 'µg/m³'
    diff_var.long_name = 'Difference: pm25_ave[0,63,:,:] - pm25'

    if 'pm25_aqm' in pm25_ds.variables:
        aqm_var = pm25_ds.variables['pm25_aqm']
    else:
        aqm_var = pm25_ds.createVariable('pm25_aqm', 'f4', ('y', 'x'))

    aqm_var[:, :] = pm25_ave
    # Add metadata
    aqm_var.units = 'µg/m³'
    aqm_var.long_name = 'pm25_ave[0,63,:,:]'

print("Difference written to pm25.nc as variable 'pm25_diff'.")

# Create a single PDF page with both histograms overlaid
# Compute shared bin edges
min_val = 0.0   #min(np.min(obs_valid), np.min(pred_valid))
max_val = 500.0 #max(np.max(obs_valid), np.max(pred_valid))
bins = np.linspace(min_val, max_val, 51)  # 50 bins

# Histogram data only (no plotting yet)
obs_hist, _ = np.histogram(obs_valid, bins=bins)
pred_hist, _ = np.histogram(pred_valid, bins=bins)

# Compute bar centers
bin_centers = 0.5 * (bins[1:] + bins[:-1])
bar_width = (bins[1] - bins[0]) * 0.4  # smaller than bin width

# Create the PDF
with PdfPages('pm25_comparison.pdf') as pdf:
    plt.figure(figsize=(10, 6))

    plt.bar(bin_centers - bar_width/2, obs_hist, width=bar_width, color='lightcoral', label='Observed (pm25)', edgecolor='k', alpha=0.8)
    plt.bar(bin_centers + bar_width/2, pred_hist, width=bar_width, color='skyblue', label='Predicted (pm25_ave)', edgecolor='k', alpha=0.8)

    plt.title('PM2.5 Comparison: Observed vs Predicted')
    plt.xlabel('PM2.5 Concentration (µg/m³)')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)

    pdf.savefig()
    plt.close()

print("Side-by-side PDF histogram saved as 'pm25_comparison.pdf'")


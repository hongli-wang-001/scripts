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
        del pm25_ds.variables['pm25_diff']

    # Define the new variable in pm25.nc
    diff_var = pm25_ds.createVariable('pm25_diff', 'f4', ('y', 'x'))
    diff_var[:, :] = full_diff 

    # Add metadata
    diff_var.units = 'µg/m³'
    diff_var.long_name = 'Difference: pm25_ave[0,63,:,:] - pm25'

print("Difference written to pm25.nc as variable 'pm25_diff'.")

# Create a single PDF page with both histograms overlaid
with PdfPages('pm25_comparison.pdf') as pdf:
    plt.figure(figsize=(8, 6))

    #plt.hist(obs_valid, bins=50, color='lightcoral', edgecolor='k', alpha=0.6, label='Observed (pm25)')
    #plt.hist(pred_valid, bins=50, color='skyblue', edgecolor='k', alpha=0.6, label='Predicted (pm25_ave[0,63,:,:])')
    # Compute shared bin edges
    min_val = min(np.min(obs_valid), np.min(pred_valid))
    max_val = 500.0 #max(np.max(obs_valid), np.max(pred_valid))
    bins = np.linspace(min_val, max_val, 51)  # 50 bins

    # Use the same bins for both histograms
    plt.hist(obs_valid, bins=bins, color='lightcoral', edgecolor='k', alpha=0.6, label='Observed (pm25)')
    plt.hist(pred_valid, bins=bins, color='skyblue', edgecolor='k', alpha=0.6, label='Predicted (pm25_ave[0,63,:,:])')

    plt.title('PM2.5 Comparison: Observed vs Predicted')
    plt.xlabel('PM2.5 Concentration (µg/m³)')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    plt.grid(True)

    pdf.savefig()
    plt.close()

print("PDF plot saved to 'pm25_comparison.pdf'")

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Open fcst.nc (read-only) and bc.nc (append mode)
with nc.Dataset('fcst.nc', 'r') as fcst_ds, nc.Dataset('bc.nc', 'a') as bc_ds:
    # Read the specific slice from fcst.nc
    eci_ave = fcst_ds.variables['aeci'][0, 63,  :, :]  # Shape: [225, 393]
    ecj_ave = fcst_ds.variables['aecj'][0, 63,  :, :]  # Shape: [225, 393]
    bc_ave = eci_ave + ecj_ave

    ak = fcst_ds.getncattr('ak')  # or ak = ds.ak
    #print("ak:", ak)
    tmp   = fcst_ds.variables['tmp'][0, 63,  :, :] 
    dpres = fcst_ds.variables['dpres'][0, :,  :, :]
    # dpres has shape (level, lat, lon)
    surface_pressure = np.sum(dpres, axis=0) + ak[0]  # sum along the 'level' dimension
    print(f"Mean surface pressure: {np.mean(surface_pressure):.3f} Pa")
    # Compute dry air density: ρ = p / (R * T)
    r_dry_air = 287.05  # J/(kg·K), specific gas constant for dry air
    dry_air_dens = surface_pressure / (tmp * r_dry_air)
    print(f"Mean dry air density: {np.mean(dry_air_dens):.3f} kg/m³")

    # ug/Kg to ug/M^3
    bc_ave = bc_ave*dry_air_dens

    # Read bc from bc.nc
    bc = bc_ds.variables['bc'][:, :]  # Shape: [225, 393]

    # Flatten and mask invalid (NaN or fill value -1)
    pred = bc_ave.flatten()
    obs = bc.flatten()

    valid_mask = (~np.isnan(pred)) & (~np.isnan(obs)) & (obs != -1) & (pred != -1)
    pred_valid = pred[valid_mask]
    obs_valid = obs[valid_mask]

    pred_mean = np.mean(pred_valid)
    obs_mean =  np.mean(obs_valid)

    # Print statistics
    print("Statistics for bc_ave[0,:,:] & bc")
    print(f"  Mean_aqm    : {pred_mean:.4f}")
    print(f"  Mean_obs    : {obs_mean:.4f}")

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
    print("Statistics for bc_diff (bc_ave[0,:,:] - bc):")
    print(f"  Valid data points: {num_valid_points}")
    print(f"  Mean      : {mean_diff:.4f}")
    print(f"  Median    : {median_diff:.4f}")
    print(f"  Std Dev   : {std_diff:.4f}")
    print(f"  RMS       : {rms_diff:.4f}")
    print(f"  IOA       : {ioa:.4f}")
    print(f"  COR       : {correlation:.4f}")

    # Save the full difference field into bc.nc
    full_diff = bc_ave - bc

    # If variable already exists, delete it to avoid conflict
    if 'bc_diff' in bc_ds.variables:
        diff_var = bc_ds.variables['bc_diff']
    else:
        diff_var = bc_ds.createVariable('bc_diff', 'f4', ('y', 'x'))

    diff_var[:, :] = full_diff 

    # Add metadata
    diff_var.units = ''
    diff_var.long_name = 'Difference: bc_aqm[0,:,:] - bc'

    if 'bc_aqm' in bc_ds.variables:
        aqm_var = bc_ds.variables['bc_aqm']
    else:
        aqm_var = bc_ds.createVariable('bc_aqm', 'f4', ('y', 'x'))

    aqm_var[:, :] = bc_ave
    # Add metadata
    aqm_var.units = 'ug/M^3'
    aqm_var.long_name = 'bc_aqm[0,:,:]'

    if 'rho_aqm' in bc_ds.variables:
        rho_var = bc_ds.variables['rho_aqm']
    else:
        rho_var = bc_ds.createVariable('rho_aqm', 'f4', ('y', 'x'))

    rho_var[:, :] = dry_air_dens
    # Add metadata
    rho_var.units = 'Kg/M^3'
    rho_var.long_name = 'dry air density[:,:]'

print("Difference written to bc.nc as variable 'bc_diff'.")

# Create a single PDF page with both histograms overlaid
# Compute shared bin edges
min_val = 0.0   #min(np.min(obs_valid), np.min(pred_valid))
max_val = 4.0 #max(np.max(obs_valid), np.max(pred_valid))
bins = np.linspace(min_val, max_val, 51)  # 50 bins

# Histogram data only (no plotting yet)
obs_hist, _ = np.histogram(obs_valid, bins=bins)
pred_hist, _ = np.histogram(pred_valid, bins=bins)

# Compute bar centers
bin_centers = 0.5 * (bins[1:] + bins[:-1])
bar_width = (bins[1] - bins[0]) * 0.4  # smaller than bin width

# Create the PDF
with PdfPages('bc_comparison.pdf') as pdf:
    plt.figure(figsize=(10, 6))

    plt.bar(bin_centers - bar_width/2, obs_hist, width=bar_width, color='lightcoral', label='AERONET AOD', edgecolor='k', alpha=0.8)
    plt.bar(bin_centers + bar_width/2, pred_hist, width=bar_width, color='skyblue', label='AQM AOD', edgecolor='k', alpha=0.8)

    plt.title('AOD Comparison: AERONET Observed vs Predicted')
    plt.xlabel('AOD')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)

    pdf.savefig()
    plt.close()

print("Side-by-side PDF histogram saved as 'bc_comparison.pdf'")


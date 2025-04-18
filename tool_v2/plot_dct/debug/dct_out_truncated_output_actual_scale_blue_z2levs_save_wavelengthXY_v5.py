import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import dct, idct
from netCDF4 import Dataset

# Load the NetCDF file and read the variable
netcdf_file = 'file1.nc'
variable_name = 'anh4k'

# Open the NetCDF file
with Dataset(netcdf_file, 'r') as ds:
    # Read the 4D variable (Time, zaxis_1, yaxis_1, xaxis_1)
    data = ds.variables[variable_name][:]
    # Assuming Time dimension is 1, extract the relevant data
    data = data[0, :, :, :]  # Shape: (zaxis_1, yaxis_1, xaxis_1)

# Parameters
z_dim, y_dim, x_dim = data.shape
z_levels = [63, 58]  # Focus on these z levels
min_wavelength = 2    # Minimum wavelength
max_wavelength = x_dim / 2  # Maximum wavelength (half of grid size)

def apply_dct_2d(array_2d):
    return dct(dct(array_2d.T, norm='ortho').T, norm='ortho')

def apply_idct_2d(array_2d_dct):
    return idct(idct(array_2d_dct.T, norm='ortho').T, norm='ortho')

def truncate_dct(dct_data, min_wavelength, max_wavelength):
    """ Truncate DCT coefficients based on specified min and max wavelengths """
    nrows, ncols = dct_data.shape
    freq_limit_min = 1 / min_wavelength
    freq_limit_max = 1 / max_wavelength
    
    # Frequency indices
    u = np.fft.fftfreq(nrows, d=1.0) * nrows
    v = np.fft.fftfreq(ncols, d=1.0) * ncols
    u_freq, v_freq = np.meshgrid(u, v, indexing='ij')
    
    # Calculate frequency magnitude
    freq_magnitude = np.sqrt(u_freq**2 + v_freq**2)
    
    # Create masks for truncation
    mask = (freq_magnitude >= freq_limit_min) & (freq_magnitude <= freq_limit_max)
    
    dct_data_truncated = np.copy(dct_data)
    dct_data_truncated[~mask] = 0
    
    return dct_data_truncated, u_freq, v_freq

def calculate_wavelengths(u_freq, v_freq):
    """ Calculate the wavelengths for each frequency component """
    freq_magnitude = np.sqrt(u_freq**2 + v_freq**2)
    
    # Avoid zero frequencies
    with np.errstate(divide='ignore', invalid='ignore'):
        wavelengths = np.where(freq_magnitude == 0, np.inf, 1 / freq_magnitude)
    
    # Filter out inf values
    valid_wavelengths = wavelengths[np.isfinite(wavelengths)]
    
    min_wavelength = np.min(valid_wavelengths) if valid_wavelengths.size > 0 else np.nan
    max_wavelength = np.max(valid_wavelengths) if valid_wavelengths.size > 0 else np.nan
    
    return min_wavelength, max_wavelength

# Plotting
plt.figure(figsize=(18, 12))

for i, z_level in enumerate(z_levels):
    # Plot the original data for the given z level
    plt.subplot(len(z_levels), 4, i * 4 + 1)
    plt.title(f'Original (z={z_level})')
    plt.imshow(data[z_level, :, :], cmap='Blues', vmin=np.min(data), vmax=np.max(data))
    plt.colorbar()

    # Apply DCT
    dct_data = apply_dct_2d(data[z_level, :, :])
    
    # Apply truncation based on new wavelength constraints
    dct_data_truncated, u_freq, v_freq = truncate_dct(dct_data, min_wavelength, max_wavelength)
    
    plt.subplot(len(z_levels), 4, i * 4 + 2)
    plt.title(f'DCT Coefficients (z={z_level})')
    plt.imshow(np.log(np.abs(dct_data) + 1e-10), cmap='Blues')  # Log scale for visibility
    plt.colorbar()

    plt.subplot(len(z_levels), 4, i * 4 + 3)
    plt.title(f'Reconstructed (z={z_level})')
    plt.imshow(apply_idct_2d(dct_data), cmap='Blues', vmin=np.min(data), vmax=np.max(data))
    plt.colorbar()

    plt.subplot(len(z_levels), 4, i * 4 + 4)
    plt.title(f'Truncated Reconstruction (z={z_level})')
    plt.imshow(apply_idct_2d(dct_data_truncated), cmap='Blues', vmin=np.min(data), vmax=np.max(data))
    plt.colorbar()

    # Calculate and print wavelengths
    min_wavelength, max_wavelength = calculate_wavelengths(u_freq, v_freq)

    print(f"\nFor z={z_level} with min wavelength {min_wavelength} and max wavelength {max_wavelength}:")
    print(f"  Min Wavelength X: {min_wavelength:.2f}")
    print(f"  Max Wavelength X: {max_wavelength:.2f}")
    print(f"  Min Wavelength Y: {min_wavelength:.2f}")
    print(f"  Max Wavelength Y: {max_wavelength:.2f}")

# Save the figure
plt.tight_layout()
plt.savefig('dct_truncated_reconstruction_actual_scale_min_max_wavelength.png')
print("DCT truncated reconstruction with actual scale and min/max wavelength figure saved as 'dct_truncated_reconstruction_actual_scale_min_max_wavelength.png'.")

# Optionally show the plots
plt.show()


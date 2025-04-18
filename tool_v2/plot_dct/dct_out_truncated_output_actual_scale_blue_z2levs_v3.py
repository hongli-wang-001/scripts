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
truncations = [10, 20, 30]  # Truncate the DCT coefficients at these levels

# Custom color map from white to blue
cmap = plt.get_cmap('Blues')
cmap.set_under('white')

# Function to apply 2D DCT and IDCT
def apply_dct_2d(array_2d):
    return dct(dct(array_2d.T, norm='ortho').T, norm='ortho')

def apply_idct_2d(array_2d_dct):
    return idct(idct(array_2d_dct.T, norm='ortho').T, norm='ortho')

def get_dct_wavelengths(dct_data):
    """ Calculate spatial scales for DCT data """
    nrows, ncols = dct_data.shape
    u = np.arange(nrows)
    v = np.arange(ncols)
    
    # Calculate wavelength based on DCT coefficients
    wavelength_x = np.zeros(ncols)
    wavelength_y = np.zeros(nrows)
    
    # Wavelength for each frequency component
    for i in range(nrows):
        for j in range(ncols):
            # Spatial frequencies
            freq_x = j if j < ncols // 2 else j - ncols
            freq_y = i if i < nrows // 2 else i - nrows
            wavelength_x[j] = np.sqrt(freq_x**2 + freq_y**2) if (freq_x**2 + freq_y**2) != 0 else np.inf
            wavelength_y[i] = np.sqrt(freq_x**2 + freq_y**2) if (freq_x**2 + freq_y**2) != 0 else np.inf

    return wavelength_x, wavelength_y

# Plotting
plt.figure(figsize=(24, 12))

for i, z_level in enumerate(z_levels):
    # Plot the original data for the given z level
    plt.subplot(len(z_levels), 5, i * 5 + 1)
    plt.title(f'Original (z={z_level})')
    plt.imshow(data[z_level, :, :], cmap=cmap, vmin=np.min(data), vmax=np.max(data))
    plt.colorbar()

    # Apply DCT
    dct_data = apply_dct_2d(data[z_level, :, :])
    
    plt.subplot(len(z_levels), 5, i * 5 + 2)
    plt.title(f'DCT Coefficients (z={z_level})')
    plt.imshow(np.log(np.abs(dct_data) + 1e-10), cmap=cmap)  # Log scale for visibility, custom cmap
    plt.colorbar()

    plt.subplot(len(z_levels), 5, i * 5 + 3)
    plt.title(f'Reconstructed (z={z_level})')
    plt.imshow(apply_idct_2d(dct_data), cmap=cmap, vmin=np.min(data), vmax=np.max(data))
    plt.colorbar()

    for trunc in truncations:
        # Truncate DCT coefficients and reconstruct
        dct_data_truncated = dct_data.copy()
        dct_data_truncated[trunc:, :] = 0
        dct_data_truncated[:, trunc:] = 0
        idct_data_truncated = apply_idct_2d(dct_data_truncated)
        
        # Calculate spatial scales
        wavelength_x, wavelength_y = get_dct_wavelengths(dct_data_truncated)
        avg_wavelength_x = np.mean(wavelength_x[wavelength_x < np.inf])
        avg_wavelength_y = np.mean(wavelength_y[wavelength_y < np.inf])
        min_wavelength_x = np.min(wavelength_x[wavelength_x < np.inf])
        max_wavelength_x = np.max(wavelength_x[wavelength_x < np.inf])
        min_wavelength_y = np.min(wavelength_y[wavelength_y < np.inf])
        max_wavelength_y = np.max(wavelength_y[wavelength_y < np.inf])
        
        plt.subplot(len(z_levels), 5, i * 5 + 4)
        plt.title(f'Trunc {trunc} (z={z_level})')
        plt.imshow(idct_data_truncated, cmap=cmap, vmin=np.min(data), vmax=np.max(data))
        plt.colorbar()
        plt.text(0.05, 0.95, f'Wavelength X: {avg_wavelength_x:.2f}', fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.7), transform=plt.gca().transAxes)
        plt.text(0.05, 0.85, f'Wavelength Y: {avg_wavelength_y:.2f}', fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.7), transform=plt.gca().transAxes)
        plt.text(0.05, 0.75, f'Min Wavelength X: {min_wavelength_x:.2f}', fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.7), transform=plt.gca().transAxes)
        plt.text(0.05, 0.65, f'Max Wavelength X: {max_wavelength_x:.2f}', fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.7), transform=plt.gca().transAxes)
        plt.text(0.05, 0.55, f'Min Wavelength Y: {min_wavelength_y:.2f}', fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.7), transform=plt.gca().transAxes)
        plt.text(0.05, 0.45, f'Max Wavelength Y: {max_wavelength_y:.2f}', fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.7), transform=plt.gca().transAxes)

        # Print wavelengths to screen
        print(f"\nFor z={z_level} and truncation={trunc}:")
        print(f"  Avg Wavelength X: {avg_wavelength_x:.2f}")
        print(f"  Avg Wavelength Y: {avg_wavelength_y:.2f}")
        print(f"  Min Wavelength X: {min_wavelength_x:.2f}")
        print(f"  Max Wavelength X: {max_wavelength_x:.2f}")
        print(f"  Min Wavelength Y: {min_wavelength_y:.2f}")
        print(f"  Max Wavelength Y: {max_wavelength_y:.2f}")

# Save the figure
plt.tight_layout()
plt.savefig('dct_truncation_reconstruction_z63_z58_with_correct_scales.png')
print("DCT truncation reconstruction for z=63 and z=58 with correct scales figure saved as 'dct_truncation_reconstruction_z63_z58_with_correct_scales.png'.")

# Optionally show the plots
plt.show()


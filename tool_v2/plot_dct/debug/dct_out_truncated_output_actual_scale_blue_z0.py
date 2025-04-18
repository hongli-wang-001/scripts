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
z_level = 0  # Show only z=0
truncations = [10, 20, 30]  # Truncate the DCT coefficients at these levels

# Custom color map from white to blue
cmap = plt.get_cmap('Blues')
cmap.set_under('white')

# Function to apply 2D DCT and IDCT
def apply_dct_2d(array_2d):
    return dct(dct(array_2d.T, norm='ortho').T, norm='ortho')

def apply_idct_2d(array_2d_dct):
    return idct(idct(array_2d_dct.T, norm='ortho').T, norm='ortho')

def get_spatial_scale(dct_data, truncation_level):
    # Spatial frequency calculation
    freq_x = np.fft.fftfreq(dct_data.shape[1])
    freq_y = np.fft.fftfreq(dct_data.shape[0])
    freq_x, freq_y = np.meshgrid(freq_x, freq_y)
    wavelength_x = 1 / (np.abs(freq_x) + 1e-10)
    wavelength_y = 1 / (np.abs(freq_y) + 1e-10)
    
    # Calculate average wavelength for the given truncation level
    avg_wavelength_x = np.mean(wavelength_x[:truncation_level, :])
    avg_wavelength_y = np.mean(wavelength_y[:, :truncation_level])
    
    return avg_wavelength_x, avg_wavelength_y

# Plotting
plt.figure(figsize=(24, 12))

# Plot the original, DCT coefficients, and reconstructed for z=0
plt.subplot(1, 4, 1)
plt.title('Original (z=0)')
plt.imshow(data[z_level, :, :], cmap=cmap, vmin=np.min(data), vmax=np.max(data))
plt.colorbar()

# Apply DCT
dct_data = apply_dct_2d(data[z_level, :, :])

plt.subplot(1, 4, 2)
plt.title('DCT Coefficients (z=0)')
plt.imshow(np.log(np.abs(dct_data) + 1e-10), cmap=cmap)  # Log scale for visibility, custom cmap
plt.colorbar()

plt.subplot(1, 4, 3)
plt.title('Reconstructed (z=0)')
plt.imshow(apply_idct_2d(dct_data), cmap=cmap, vmin=np.min(data), vmax=np.max(data))
plt.colorbar()

# Truncate DCT coefficients and reconstruct for different truncations
for trunc in truncations:
    dct_data_truncated = dct_data.copy()
    dct_data_truncated[trunc:, :] = 0
    dct_data_truncated[:, trunc:] = 0
    idct_data_truncated = apply_idct_2d(dct_data_truncated)
    
    # Calculate spatial scales
    avg_wavelength_x, avg_wavelength_y = get_spatial_scale(dct_data_truncated, trunc)
    
    plt.subplot(1, 4, 4)
    plt.title(f'Trunc {trunc} (z=0)')
    plt.imshow(idct_data_truncated, cmap=cmap, vmin=np.min(data), vmax=np.max(data))
    plt.colorbar()
    plt.text(0.05, 0.95, f'Wavelength X: {avg_wavelength_x:.2f}', fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.7), transform=plt.gca().transAxes)
    plt.text(0.05, 0.85, f'Wavelength Y: {avg_wavelength_y:.2f}', fontsize=10, color='red', bbox=dict(facecolor='white', alpha=0.7), transform=plt.gca().transAxes)

# Save the figure
plt.tight_layout()
plt.savefig('dct_truncation_reconstruction_with_scales_all_plots.png')
print("DCT truncation reconstruction with updated color maps and scales figure saved as 'dct_truncation_reconstruction_with_scales_all_plots.png'.")

# Optionally show the plots
plt.show()


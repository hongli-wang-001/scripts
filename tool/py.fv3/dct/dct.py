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
scales_to_plot = [0, z_dim // 2, z_dim - 1]  # Example: first, middle, last z levels

# Function to perform 2D DCT and IDCT
def apply_dct_2d(array_2d):
    return dct(dct(array_2d.T, norm='ortho').T, norm='ortho')

def apply_idct_2d(array_2d_dct):
    return idct(idct(array_2d_dct.T, norm='ortho').T, norm='ortho')

# Plotting
plt.figure(figsize=(15, 10))

# Process each z dimension level
for idx in scales_to_plot:
    plt.subplot(3, 3, scales_to_plot.index(idx) * 3 + 1)
    plt.title(f'Original (z={idx})')
    plt.imshow(data[idx, :, :], cmap='gray')
    plt.colorbar()
    
    # Apply DCT
    dct_data = apply_dct_2d(data[idx, :, :])
    
    # Apply IDCT (for reconstruction)
    idct_data = apply_idct_2d(dct_data)
    
    plt.subplot(3, 3, scales_to_plot.index(idx) * 3 + 2)
    plt.title(f'DCT Coefficients (z={idx})')
    plt.imshow(np.log(np.abs(dct_data) + 1e-10), cmap='gray')  # Log scale for visibility
    plt.colorbar()
    
    plt.subplot(3, 3, scales_to_plot.index(idx) * 3 + 3)
    plt.title(f'Reconstructed (z={idx})')
    plt.imshow(idct_data, cmap='gray')
    plt.colorbar()

# Save the figure
plt.tight_layout()
plt.savefig('dct_reconstruction.png')
print("DCT reconstruction figure saved as 'dct_reconstruction.png'.")

# Optionally show the plots
plt.show()


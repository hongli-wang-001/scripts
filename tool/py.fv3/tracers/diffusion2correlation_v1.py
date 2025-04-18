import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.signal import correlate2d

# Generate synthetic 2D spatial data
np.random.seed(0)
data = np.random.rand(100, 100) * 100  # Random 2D data representing spatial field

# Apply Gaussian filter to model diffusion
sigma = 5.0  # Standard deviation for Gaussian kernel
smoothed_data = gaussian_filter(data, sigma=sigma)

# Compute correlation matrix of the smoothed data
correlation = correlate2d(smoothed_data, smoothed_data, mode='full')

# Plot original and smoothed data
fig, ax = plt.subplots(1, 3, figsize=(18, 6))

# Plot original data
c1 = ax[0].imshow(data, cmap='viridis', origin='lower')
ax[0].set_title('Original Data')
fig.colorbar(c1, ax=ax[0])

# Plot smoothed data
c2 = ax[1].imshow(smoothed_data, cmap='viridis', origin='lower')
ax[1].set_title('Smoothed Data (Diffusion)')
fig.colorbar(c2, ax=ax[1])

# Plot the correlation matrix
c3 = ax[2].imshow(correlation, cmap='viridis', origin='lower')
ax[2].set_title('Correlation Matrix')
fig.colorbar(c3, ax=ax[2])

# Save the figure
plt.savefig('diffusion_correlation_analysis.png', dpi=300)

# Show the figure
plt.show()


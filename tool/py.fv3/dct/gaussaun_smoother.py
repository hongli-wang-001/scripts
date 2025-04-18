import numpy as np
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt

# Create a sample 2D array
data = np.random.rand(100, 100) * 100

# Apply Gaussian filter
smoothed_data = gaussian_filter(data, sigma=2)

# Plot original and smoothed data
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.title("Original Data")
plt.imshow(data, cmap='viridis')
plt.colorbar()

plt.subplot(1, 2, 2)
plt.title("Smoothed Data")
plt.imshow(smoothed_data, cmap='viridis')
plt.colorbar()

# Save the figure
plt.savefig('smoothed_data_gaussian_filter.png', dpi=300)

plt.show()


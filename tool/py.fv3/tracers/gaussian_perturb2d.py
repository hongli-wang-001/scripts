import numpy as np
import scipy.ndimage as ndimage
import matplotlib.pyplot as plt

def apply_2d_gaussian_filter(original_field, sigma):
    """
    Apply a 2D Gaussian filter to a 2D field with different sigmas for x and y directions.

    Parameters:
    - original_field: 2D numpy array of the original field
    - sigma: Tuple of standard deviations (sigma_x, sigma_y)

    Returns:
    - 2D numpy array with applied 2D Gaussian filter
    """
    filtered_field = ndimage.gaussian_filter(original_field, sigma=sigma)
    return filtered_field

# Example 2D field (replace with your actual data)
original_field = np.random.rand(100, 100)

# Define 2D sigma (standard deviations for x and y directions)
sigma = (2.0, 5.0)

# Apply the 2D Gaussian filter
filtered_field = apply_2d_gaussian_filter(original_field, sigma)

# Plot the original and filtered fields
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
axes[0].imshow(original_field, cmap='viridis')
axes[0].set_title('Original Field')
axes[1].imshow(filtered_field, cmap='viridis')
axes[1].set_title('Filtered Field')
plt.tight_layout()
plt.show()


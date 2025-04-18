import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

# Example data
data = np.random.random((100, 100))

# Define sigma values to test
sigma_values = [1, 3, 5]

plt.figure(figsize=(12, 4))

for i, sigma in enumerate(sigma_values):
    filtered_data = gaussian_filter(data, sigma=sigma)
    plt.subplot(1, len(sigma_values), i + 1)
    plt.imshow(filtered_data, cmap='viridis')
    plt.title(f'Sigma = {sigma}')
    plt.axis('off')

plt.tight_layout()
plt.show()


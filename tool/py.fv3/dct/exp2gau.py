import numpy as np
import matplotlib.pyplot as plt

# Generate sample data with a positive skew (e.g., exponential distribution)
np.random.seed(0)
data = np.random.exponential(scale=2, size=1000)

# Apply Square Root Transformation
transformed_data = np.sqrt(data)

# Create the figure and axes
fig, ax = plt.subplots(1, 2, figsize=(12, 6))

# Plot original data
ax[0].hist(data, bins=30, alpha=0.7, color='blue', edgecolor='black')
ax[0].set_title('Original Data')
ax[0].set_xlabel('Value')
ax[0].set_ylabel('Frequency')

# Plot transformed data
ax[1].hist(transformed_data, bins=30, alpha=0.7, color='green', edgecolor='black')
ax[1].set_title('Square Root Transformed Data')
ax[1].set_xlabel('Transformed Value')
ax[1].set_ylabel('Frequency')

# Adjust layout to avoid overlapping
plt.tight_layout()

# Save the figure to a file
plt.savefig('square_root_transformation.png')

# Optionally, you can display the figure if needed
# plt.show()


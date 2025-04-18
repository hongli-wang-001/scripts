import numpy as np
import matplotlib.pyplot as plt

# Example data
x_data = np.random.randn(1000)
y_data = np.random.randn(1000)

# Plot with hexbin
plt.figure(figsize=(10, 6))
hexbin = plt.hexbin(x_data, y_data, gridsize=50, cmap='Blues', mincnt=1)

# Custom scaling (e.g., top 90% density)
max_density = np.percentile(hexbin.get_array(), 90)
plt.colorbar(label='Density', extend='max')

# Customize the color limits
plt.clim(0, max_density)  # Set color limits based on the percentile

plt.xlabel('X Data')
plt.ylabel('Y Data')
plt.title('Hexbin Plot with Custom Density Scaling')
plt.grid(True)

# Save the figure
plt.savefig('hexbin_plot_with_density_scaling.png', dpi=300)  # Save with high resolution

# Optionally, you can also display the plot if needed
# plt.show()

print('Figure saved as hexbin_plot_with_density_scaling.png')


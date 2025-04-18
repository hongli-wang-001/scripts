import numpy as np
import matplotlib.pyplot as plt

# Define the parameters of the normal distribution
mu = 0  # mean
sigma = 1  # standard deviation

# Generate random numbers from the normal distribution
x = np.random.normal(mu, sigma, 10000)

# Plot the PDF of the normal distribution
plt.hist(x, bins=50, density=True, alpha=0.6, color='g')
plt.xlabel('Value')
plt.ylabel('Probability density')
plt.title('Normal distribution PDF')
plt.show()

import numpy as np
import matplotlib.pyplot as plt

# Load the data from the files
exp_data = np.loadtxt("exponential_data.dat")
gauss_data = np.loadtxt("gaussian_data.dat")

# Plot the exponential data
plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.hist(exp_data, bins=30, density=True, alpha=0.7, color='blue', label='Exponential Data')
plt.title('Exponential Distribution')
plt.xlabel('Value')
plt.ylabel('Density')
plt.legend()

# Plot the Gaussian data
plt.subplot(1, 2, 2)
plt.hist(gauss_data, bins=30, density=True, alpha=0.7, color='green', label='Gaussian Data')
plt.title('Gaussian Distribution')
plt.xlabel('Value')
plt.ylabel('Density')
plt.legend()

# Show the plots
plt.tight_layout()
plt.savefig("distributions_plot.png")
plt.show()

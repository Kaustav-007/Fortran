import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps  # Correct import

# Load data from the .dat file
data1 = np.loadtxt('dist_rand_walk_1k_bin_10.dat')

# Extract x and y from the data
x = data1[:, 0]
y = data1[:, 1]

# Plot data from the file
plt.scatter(x, y, label='Bin size=10', color='blue')

# Add labels and title
plt.xlabel('Sum of random numbers')
plt.ylabel('Distribution')
plt.title('Distribution of sum of 100k random walk steps for different bin sizes normalized Gaussian')
plt.grid()

# Calculate area under the curve using Simpson's rule
area = simps(y, x)
print(f"Area under the curve: {area}")

# Show the legend
plt.legend()

# Display the plot
plt.show()


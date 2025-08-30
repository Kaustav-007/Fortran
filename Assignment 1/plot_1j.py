import numpy as np
import matplotlib.pyplot as plt

# Load data from the three .dat files
# Adjust delimiter as needed (e.g., space, tab, comma)
data1 = np.loadtxt('dist_rand_walk_bin_2.dat')
data2 = np.loadtxt('dist_rand_walk_bin_5.dat')
data3 = np.loadtxt('dist_rand_walk_bin_10.dat')

# Example: assuming each file has two columns (x, y)
# If the structure differs, adjust the slicing accordingly

# Plot data from file 1
plt.scatter(data1[:, 0], data1[:, 1], label='Bin size=2', color='blue')

# Plot data from file 2
plt.scatter(data2[:, 0], data2[:, 1], label='Bin size=5', color='green')

# Plot data from file 3
plt.scatter(data3[:, 0], data3[:, 1], label='Bin size=10', color='red')

# Add labels and title
plt.xlabel('Sum of random numbers')
plt.ylabel('Distribution')
plt.title('Distribution of sum of 10k random walk steps for diiferent bin sizes')
plt.grid()

# Show the legend
plt.legend()

# Display the plot
plt.show()


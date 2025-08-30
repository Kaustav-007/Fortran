import numpy as np
import matplotlib.pyplot as plt

# Load data from the three .dat files
# Adjust delimiter as needed (e.g., space, tab, comma)
data1 = np.loadtxt('distribution.dat')
data2 = np.loadtxt('distribution1.dat')
data3 = np.loadtxt('distribution2.dat')

# Example: assuming each file has two columns (x, y)
# If the structure differs, adjust the slicing accordingly

# Plot data from file 1
plt.scatter(data1[:, 0], data1[:, 1], label='Bin size=2', color='blue')

# Plot data from file 2
plt.scatter(data2[:, 0], data2[:, 1], label='Bin size=1', color='green')

# Plot data from file 3
plt.scatter(data3[:, 0], data3[:, 1], label='Bin size=0.5', color='red')

# Add labels and title
plt.xlabel('Sum of random numbers')
plt.ylabel('Distribution')
plt.title('Distribution of sum of 10k random numbers')
plt.grid()

# Show the legend
plt.legend()

# Display the plot
plt.show()


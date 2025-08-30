import numpy as np
import matplotlib.pyplot as plt

# Load data from the three .dat files
# Adjust delimiter as needed (e.g., space, tab, comma)
data1 = np.loadtxt('dist_range(-1 TO +1)')
data2 = np.loadtxt('dist_range_100k(-1 TO +1)')


# Example: assuming each file has two columns (x, y)
# If the structure differs, adjust the slicing accordingly

# Plot data from file 1
plt.scatter(data1[:, 0], data1[:, 1], label='10k times', color='blue')

# Plot data from file 2
plt.scatter(data2[:, 0], data2[:, 1], label='100k times', color='green')


# Add labels and title
plt.xlabel('Sum of random numbers')
plt.ylabel('Distribution')
plt.title('Distribution of sum of 10k and 100k random numbers')
plt.grid()
# Show the legend
plt.legend()

# Display the plot
plt.show()


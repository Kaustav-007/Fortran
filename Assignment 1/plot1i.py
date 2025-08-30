import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps  # Correct import

# Load data from the three .dat files
# Adjust delimiter as needed (e.g., space, tab, comma)
data1 = np.loadtxt("dist_new.dat")
#data2= np.loadtxt("dist_rand_walk_bin_0.5.dat")

x=data1[:,0]
y=data1[:,1]

# Example: assuming each file has two columns (x, y)
# If the structure differs, adjust the slicing accordingly

# Plot data from file 1
plt.scatter(x, y, label='Bin size=1', color='blue')

# Plot data from file 2
#plt.scatter(data2[:, 0], data2[:, 1], label='Bin size=0.5', color='green')


# Add labels and title
plt.xlabel('Sum of random numbers')
plt.ylabel('Distribution')
plt.title('Distribution of sum of random walk 10k step')
area = simps(y, x)
print(f"Area under the curve: {area}")

plt.grid()

# Show the legend
plt.legend()

# Display the plot
plt.show()


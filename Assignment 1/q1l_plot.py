import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import simps

# Load data from the .dat file
data1 = np.loadtxt('dist_l.dat')

# Extract x and y from the data
x = data1[:, 0]
y = data1[:, 1]
x0=0

# Define the Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

# Initial guess for the parameters [amplitude, mean, standard deviation]
initial_guess = [max(y), np.mean(x), np.std(x)]

# Fit the Gaussian function to the data
popt, pcov = curve_fit(gaussian, x, y, p0=initial_guess)

# Plot data and the fitted curve
plt.scatter(x, y, label='Data', color='blue')
plt.plot(x, gaussian(x, *popt), label='Fitted Gaussian', color='red')

# Add labels and title
plt.xlabel('Sum of random numbers')
plt.ylabel('Distribution')
plt.title('Distribution of sum of 100k random walk steps for bin size= 10 normalized Gaussian')
plt.grid()

# Calculate area under the curve using Simpson's rule
area = simps(y, x)
print(f"Area under the curve: {area}")

# Show the legend
plt.legend()

# Display the plot
plt.show()

# Print the fitted parameters
print(f"Fitted parameters: Amplitude = {popt[0]},  Std Dev = {popt[2]}")


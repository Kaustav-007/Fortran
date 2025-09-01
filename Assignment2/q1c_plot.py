import numpy as np
import matplotlib.pyplot as plt

# Load data from the Fortran output file
data = np.loadtxt('trap_sin.dat')
n = data[:, 0]  # Number of intervals
errors = data[:, 1]  # Errors

# Log-log plot
plt.loglog(n, errors, marker='o', linestyle='-', label='Error')
#plt.xlim(100, 100000)  # x-range from -15 t
#plt.ylim(1e-11, 1e-2)

# Fit a line to log-log data

# Display the fitted line
plt.loglog(n, errors, linestyle='--', )

plt.xlabel('Number of Intervals (n)')
plt.ylabel('Error')
plt.title('Log-Log Plot of Error vs. Number of Intervals')
plt.legend
plt.grid()

# Show the plot
plt.show()


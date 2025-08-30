import numpy as np
import matplotlib.pyplot as plt
import scienceplots
plt.style.use(['science', 'grid'])



# Load data from files
data1 = np.loadtxt('euler1.dat')
data2 = np.loadtxt('euler_mod1.dat')
data3 = np.loadtxt('euler_mod2.dat')
data4=  np.loadtxt('RK4.dat')

# Extract x and y values
x1, y1 = data1[:, 0], data1[:, 1]
x2, y2 = data2[:, 0], data2[:, 1]
x3, y3 = data3[:, 0], data3[:, 1]
x4, y4 = data4[:, 0], data4[:, 1]

# Define x values for y = tan(x), avoiding asymptotes
x_exact = np.linspace(0, 1.55)  # Adjust range to avoid singularities
y_exact = np.tan(x_exact)

# Plot the data
plt.figure(figsize=(10, 6))
plt.plot(x1, y1, 'r', label='Euler Method')
plt.plot(x2, y2, '#FFFF00',lw=3, label='Modified Euler Method')
plt.plot(x3, y3, 'b', label='Improved Euler Method')
plt.plot(x4, y4, 'g', label='RK4 Method')
plt.plot(x_exact, y_exact, 'k', label='y = tan(x)')  # Add tan(x) plot

# Formatting the plot
plt.xlabel('x')
plt.ylabel('y')
plt.title('Comparison of Euler Methods with y = tan(x)',fontsize=30)
plt.legend()
plt.grid(True)
#plt.ylim(-10, 10)  # Limit y-axis to avoid extreme values from tan(x)

# Show the plot
plt.show()


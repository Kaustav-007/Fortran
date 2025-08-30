import matplotlib.pyplot as plt
import numpy as np
import scienceplots
plt.style.use(['science'])

# Load data
data1 = np.loadtxt("q5_SHO.dat")
data2 = np.loadtxt("q6_SHO.dat")
data3 = np.loadtxt("q7_SHO.dat")

# Extract x and v from columns (assuming column 2 is x and column 3 is v, i.e., 1 and 2 in Python index)
x1, v1 = data1[:,1], data1[:,2]
x2, v2 = data2[:,1], data2[:,2]
x3, v3 = data3[:,1], data3[:,2]

# Plot
plt.figure(figsize=(10, 6))
plt.plot(x1, v1, label="x(0)=0.1, v(0)=1.900")
plt.plot(x2, v2, label="x(0)=0.0, v(0)=1.999")
plt.plot(x3, v3, label="x(0)=0.0, v(0)=2.00001")

# Axis labels and title
plt.xlabel("x", fontsize=14)
plt.ylabel("v", fontsize=14)
plt.title("Phase-Space : Simple Pendulum For Different Initial Conditions", fontsize=16)

# Axis limits
#plt.xlim([-3.5, 10])
#plt.ylim([-2.5, 2.5])

# Legend
plt.legend(loc='lower right', fontsize=12)

# Grid for better visibility
plt.grid()

# Show plot
#plt.tight_layout()
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

data1 = np.loadtxt('q6a.dat')
data2 = np.loadtxt('q6b.dat')
data3 = np.loadtxt('q6c.dat')

x1 = data1[:,0]
y1 = data1[:,1]
x2 = data2[:,0]
y2 = data2[:,1]
x3 = data3[:,0]
y3 = data3[:,1]

plt.plot(x1, y1,color = "black", label='Lattice Size in one Dimension = 8')
plt.plot(x2, y2,color = "blue", label='Lattice Size in one Dimension = 9')
plt.plot(x3, y3,color = "green", label='Lattice Size in one Dimension = 10')
plt.xlabel('Iteration', fontsize = 22)
plt.ylabel(r'$E/N$ Fluctuation', fontsize = 22)
plt.title(r'Energy Fluctuation vs Iteration Plot', fontsize = 24)
plt.legend(fontsize = 20, frameon = True)
plt.grid()
plt.show()

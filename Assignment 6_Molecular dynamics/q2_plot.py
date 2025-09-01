import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

data1 = np.loadtxt('energy_q2.dat')


x1 = data1[:,0]
y1 = data1[:,1]
y2 = data1[:,2]
y3 = data1[:,3]

plt.plot(x1, y1,color = "black", label='PE')
plt.plot(x1, y2,color = "blue", label='KE')
plt.plot(x1, y3,color = "green", label='Total Energy')
plt.xlabel('No of Iterations', fontsize = 22)
plt.ylabel(r'Energies', fontsize = 22)
plt.title(r'Energy vs Iteration without thermostat', fontsize = 24)
plt.legend(fontsize = 20, frameon = True)
plt.grid()
plt.show()

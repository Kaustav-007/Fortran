import numpy as np

import matplotlib.pyplot as plt

data_1=np.loadtxt('q5_ising.dat')

plt.plot(data_1[:,0],data_1[:,1],label='Energy_fluctuations',color='red')
plt.plot(data_1[:,0],data_1[:,2],label='Magnetic_fluctuations',color='green')
plt.legend()
plt.grid()

plt.xlabel('iteration')

plt.ylabel('fluctuations')


plt.show()

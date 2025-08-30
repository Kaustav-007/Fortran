import numpy as np

import matplotlib.pyplot as plt

data_1=np.loadtxt('q4_ising.dat')

plt.plot(data_1[:,0],data_1[:,1],label='Energy_fluctuations',color='red')
plt.grid()

plt.xlabel('iteration')

plt.ylabel('Energy_fluctuations')

plt.show()

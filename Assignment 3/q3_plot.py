import numpy as np

import matplotlib.pyplot as plt

data_1=np.loadtxt('q3_ising.dat')

plt.plot(data_1[:,0],data_1[:,1],label='magnetic_fluctuations',color='red')
plt.grid()

plt.xlabel('iteration')

plt.ylabel('Magnetic_fluctuations')

plt.show()

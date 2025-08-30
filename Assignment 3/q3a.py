import numpy as np

import matplotlib.pyplot as plt
plt.style.use('ggplot')       # Red and gray grid-style
#plt.style.use('seaborn')      # Seaborn-like appearance
#plt.style.use('fivethirtyeight')  # Inspired by FiveThirtyEight


data_1=np.loadtxt('q3a_ising.dat')

plt.plot(data_1[:,0],data_1[:,1],label='magnetic_fluctuations',color='red')
plt.grid()

plt.xlabel('iteration')

plt.ylabel('Energy_fluctuations')

plt.show()

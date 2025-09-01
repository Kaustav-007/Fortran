import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

# Load data from the Fortran output file
data = np.loadtxt('q4_paircorrelation.dat')
x = data[:,0]  # Number of intervals
y = data[:,1]  # Errors


plt.plot(x,y)
#plt.xlim(100, 100000)  # x-range from -15 t
#plt.ylim(1e-11, 1e-2)


plt.xlabel('r')
plt.ylabel('g(r)')
plt.title('Pair correlation plot for 3600 particles',fontsize=24)
plt.legend
plt.grid()

# Show the plot
plt.show()


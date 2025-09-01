import numpy as np
import matplotlib.pyplot as plt
import scienceplots

# Set the style for the plots
plt.style.use(['science'])
# ---------------------=================---------------------
#---------------------=================---------------------

#---------------------=================---------------------
# ---------------------=================---------------------
data = np.loadtxt('momentum_q3.dat')
data1=np.loadtxt("mom_xyz.dat")
# Extract the columns from the data
t = data[:, 0]  # Number of Iteration

p = data[:, 1]  # Total Momentum
y1=data1[:,0]
y2=data1[:,1]
y3=data1[:,2]
#total momentum plot
plt.plot(t, p, label='Total Momentum', color='blue')
plt.plot(t, y1, label=r'$p_x$', color='red')
plt.plot(t, y2, label=r'$p_y$', color='green')
plt.plot(t, y3, label=r'$p_z$', color='black')


plt.grid()
plt.title(' Momentum vs Number of Iteration ',fontsize=24)
plt.xlabel('Number of Iteration ',fontsize=22)
plt.ylabel(' Momentum',fontsize=22)
#np.polyfit(t, p, 1) #polyfit (x_array, y_array, degree)
#plt.plot(t, np.polyval(np.polyfit(t, T, 1), t), 'k-',linewidth = 2, label='Linear fit')
plt.legend()
#plt.ylim(-2.5, 2.5)
plt.show()
# ---------------------=================---------------------
# ---------------------=================---------------------

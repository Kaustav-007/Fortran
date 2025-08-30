import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')

data1 = np.loadtxt('q7_L7.dat')
data2 = np.loadtxt('q7_L8.dat')
data3 = np.loadtxt('q7_L9.dat')

x1 = data1[:,0]
y1 = data1[:,1]
x2 = data2[:,0]
y2 = data2[:,1]
x3 = data3[:,0]
y3 = data3[:,1]

plt.plot(x1, y1,'o-',color = "black", label='L=7')
plt.plot(x2, y2,'o-',color = "blue", label='L= 8')
plt.plot(x3, y3,'o-',color = "green", label='L= 9')
plt.xlabel('Temperature', fontsize = 22)
plt.ylabel('Magnetisation per Spin', fontsize = 22)
plt.title(r'Magnetization vs Temperature Plot', fontsize = 24)
plt.legend(fontsize = 20, frameon = True)
plt.grid()
plt.show()

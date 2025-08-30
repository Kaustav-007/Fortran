import numpy as np
import matplotlib.pyplot as plt
import scienceplots

plt.style.use('science')
# Load the data
data1 = np.loadtxt("q5_SHO.dat")
data2 = np.loadtxt("q6_SHO.dat")
data3 = np.loadtxt("q7_SHO.dat")
data4 = np.loadtxt("q6_v.dat")
# ---------------------=================---------------------
#---------------------=================---------------------

#---------------------=================---------------------
# ---------------------=================---------------------
#data files
#Extract position (x) and velocity (v) with respect to time
x1, v1 = data1[:,1], data1[:,2]
t1 = data1[:, 0]  # First column is time
 

x2, v2 = data2[:,1], data2[:,2]
t2 = data2[:, 0]  # First column is time


x3, v3 = data3[:,1], data3[:,2]
t3 = data3[:, 0]  # First column is time


x4, v4 = data4[:,1], data4[:,2]
t4 = data4[:, 0]  # First column is time


E1 = data1[:, 3]  # Fourth column is E
E2 = data2[:, 3]  # Fourth column is E
E3 = data3[:, 3]  # Fourth column is E
E4 = data4[:, 3]  # Fourth column is E
# ---------------------=================---------------------
# ---------------------=================---------------------


# ---------------------=================---------------------
#V vs T Plot
#plt.subplot(2,1,2)
plt.plot(t1, v1,'-', linewidth=1, label=r"$v_0 = 1.9$") 
plt.plot(t2, v2, '.-', ms = 1,linewidth=0.5, label=r"$v_0 = 1.999.$")
plt.plot(t3, v3, '--', linewidth=1.5, label=r"$v_0 = 2.00001$")
plt.plot(t4, v4, '--', linewidth=1.5, label=r"$v_0 = 2.0$")
plt.legend()
plt.xlabel(r"Time ($t$)")
plt.ylabel(r"Velocity ($v(t)$)")
plt.title("$v$ vs $t$ graph of Simple Pendulum")

plt.grid(True)
# ---------------------=================---------------------

plt.show()

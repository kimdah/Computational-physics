import numpy as np
import matplotlib.pyplot as plt

#read out w_v and N for each f as w_1,w_2,w_3,N_1,N_2,N_3

# f = 0.1
data1 = np.loadtxt('./Results/problem10_f_0.1.txt', skiprows=1)
w1 = np.array(data1[:,0])
f1 = np.array(data1[:,1])
plt.plot(w1, f1, label ="f = 0.1")

# f = 0.4
data2 = np.loadtxt('./Results/problem10_f_0.4.txt', skiprows=1)
w2= np.array(data2[:,0])
f2 = np.array(data2[:,1])
plt.plot(w2,f2, label ="f= 0.4")

# f= 0.7
data3 = np.loadtxt('./Results/problem10_f_0.7.txt', skiprows=1)
w3 = np.array(data3[:,0])
f3 = np.array(data3[:,1])
plt.plot(w3,f3, label ="f= 0.7")

plt.xlabel("Angular frequency $\omega_V$")
plt.ylabel("Fraction of escaped particles")
plt.legend()
plt.grid()
plt.savefig('fraction_vs_angfreq.pdf')
plt.show()

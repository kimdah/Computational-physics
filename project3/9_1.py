import numpy as np
import matplotlib.pyplot as plt

# Plot of z-trajectory of 1 particle

data = np.loadtxt('./Results/i_10000_d_100_p_1_pi_1_axis_z.txt', skiprows=1)
#data = data[1,:]
print(data[0,:])
t = np.array(data[:,0])
z = np.array(data[:,1])


plt.plot(t, z)
plt.title('Particle motion in z-direction for T = 100 microsec.')
plt.xlabel('time [t] = microsec.')
plt.ylabel('Position in z-direction, z(t), [z] = micrometer')
plt.savefig('9_1plot.pdf')
plt.grid()
plt.show()

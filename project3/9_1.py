import numpy as np
import matplotlib as plt


# Plot of z-trajectory of 1 particle

data = np.loadtxt()


plt.plot(t, z)
plt.title('Particle motion in z-direction for T = 100 microsec.')
plt.xlabel('time [t] = microsec.')
plt.ylabel('Position in x-direction, z(t)')
plt.savefig('9_1plot.pdf')
plt.grid()
plt.show()

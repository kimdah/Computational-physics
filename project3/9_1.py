import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Plot of z-trajectory of 1 particle

data = np.loadtxt('./Results/RK4_i_10000_d_100_p_1_pi_1_outputs_tz_pert_0_rs_0_f_0.0_w_v_0.0.txt', skiprows=1)
t = np.array(data[:,0])
z = np.array(data[:,1])

# Tick frequency on axes
ax = plt.axes()
ax.xaxis.set_major_locator(ticker.MultipleLocator(5)) # 5 = tick increment

plt.plot(t, z)

plt.title('Particle motion in z-direction for T = 100 microsec.')
plt.xlabel('time, [t] = microsec.')
plt.ylabel('Position in z-direction, z(t), [z] = micrometer')
plt.grid()



plt.savefig('./Figures/zt.pdf')
plt.show()

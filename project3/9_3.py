# Phase space plots - position vs velocity in each direction for two particles!!!
import numpy as np
import matplotlib.pyplot as plt

# WITHOUT INTERACTIONS
data0x = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_0_outputs_xv.txt', skiprows = 1)
data0y = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_0_outputs_yv.txt', skiprows = 1)
data0z = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_0_outputs_zv.txt', skiprows = 1)
# Positions:
x0 = np.array(data0x[:,0])
y0 = np.array(data0y[:,0])
z0 = np.array(data0z[:,0])
# Velocities
vx0 = np.array(data0x[:,1])
vy0 = np.array(data0y[:,1])
vz0 = np.array(data0z[:,1])

# Plot:
plt.plot(x0,v0, label='x-v_x plot')
plt.show()
plt.plot(y0,v0, label='y-v_y plot')
plt.show()
plt.plot(z0,v0, label='z-v_z plot')
plt.show()

# WITH INTERACTIONS
data1x = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_1_outputs_xv.txt', skiprows = 1)
data1y = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_1_outputs_yv.txt', skiprows = 1)
data1z = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_1_outputs_zv.txt', skiprows = 1)

# Positions
x1 = np.array(data1x[:,0])
y1 = np.array(data1y[:,0])
z1 = np.array(data1z[:,0])
# Velocities
vx1 = np.array(data1x[:,1])
vy1 = np.array(data1y[:,1])
vz1 = np.array(data1z[:,1])

# Plot
plt.plot(x1,v1, label='x-v_x w/interaction')
plt.show()
plt.plot(y1,v1, label='y-v_y w/interaction')
plt.show()
plt.plot(z1,v1, title='z-v_z w/interaction')
plt.show()

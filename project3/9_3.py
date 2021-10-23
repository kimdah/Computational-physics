# Phase space plots - position vs velocity in each direction for two particles!!!
import numpy as np
import matplotlib.pyplot as plt

# $axis $interaction $particle

# WITHOUT INTERACTIONS
data0x = np.loadtxt('./Results/RK4_i_10000_d_100_p_2_pi_0_outputs_xv_pert_0_rs_0_f_0.0_w_v_0.0.txt', skiprows = 1)
data0y = np.loadtxt('./Results/RK4_i_10000_d_100_p_2_pi_0_outputs_yv_pert_0_rs_0_f_0.0_w_v_0.0.txt', skiprows = 1)
data0z = np.loadtxt('./Results/RK4_i_10000_d_100_p_2_pi_0_outputs_zv_pert_0_rs_0_f_0.0_w_v_0.0.txt', skiprows = 1)

# ------Particle 1--------
# positions
x01 = np.array(data0x[:,0])
y01 = np.array(data0y[:,0])
z01 = np.array(data0z[:,0])
# velocities
vx01 = np.array(data0x[:,1])
vy01 = np.array(data0y[:,1])
vz01 = np.array(data0z[:,1])

# ------- Particle 2 ---------
# positions
x02 = np.array(data0x[:,2])
y02 = np.array(data0y[:,2])
z02 = np.array(data0z[:,2])
# velocities
vx02 = np.array(data0x[:,3])
vy02 = np.array(data0y[:,3])
vz02 = np.array(data0z[:,3])

# WITH INTERACTIONS

data1x = np.loadtxt('./Results/RK4_i_10000_d_100_p_2_pi_1_outputs_xv_pert_0_rs_0_f_0.0_w_v_0.0.txt', skiprows = 1)
data1y = np.loadtxt('./Results/RK4_i_10000_d_100_p_2_pi_1_outputs_yv_pert_0_rs_0_f_0.0_w_v_0.0.txt', skiprows = 1)
data1z = np.loadtxt('./Results/RK4_i_10000_d_100_p_2_pi_1_outputs_zv_pert_0_rs_0_f_0.0_w_v_0.0.txt', skiprows = 1)

# -------Particle 1 --------
# Positions
x11 = np.array(data1x[:,0])
y11 = np.array(data1y[:,0])
z11 = np.array(data1z[:,0])
# Velocities
vx11 = np.array(data1x[:,1])
vy11 = np.array(data1y[:,1])
vz11 = np.array(data1z[:,1])

# -------Particle 2 --------
# Positions
x12 = np.array(data1x[:,2])
y12 = np.array(data1y[:,2])
z12 = np.array(data1z[:,2])
# Velocities
vx12 = np.array(data1x[:,3])
vy12 = np.array(data1y[:,3])
vz12 = np.array(data1z[:,3])


# ------ PLOT---------

figx = plt.figure()
plt.title('x-v_x plot')
plt.plot(x01, vx01, label='p1 w/o ')
plt.plot(x02, vx02, label='p2 w/o')
plt.plot(x11, vx11, label='p1 w/')
plt.plot(x12, vx12, label='p2 w/')
plt.legend()
plt.savefig('./Figures/xvx.pdf')
plt.show()

figy = plt.figure()
plt.title('y-v_y plot')
plt.plot(y01, vy01, label='p1 w/o ')
plt.plot(y02, vy02, label='p2 w/o')
plt.plot(y11, vy11, label='p1 w/')
plt.plot(y12, vy12, label='p2 w/')
plt.legend()
plt.savefig('./Figures/yvy.pdf')
plt.show()

figz = plt.figure()
plt.title('z-v_z plot')
plt.plot(z01, vz01, label='p1 w/o ')
plt.plot(z02, vz02, label='p2 w/o')
plt.plot(z11, vz11, label='p1 w/')
plt.plot(z12, vz12, label='p2 w/')
plt.legend()
plt.savefig('./Figures/zvz.pdf')
plt.show()

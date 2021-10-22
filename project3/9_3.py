# Phase space plots - position vs velocity in each direction
import numpy as np
import matplotlib.pyplot as plt

# without interactions
# Positions:
data0x = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_0_outputs_xv.txt', skiprows = 1)
data0y = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_0_outputs_yv.txt', skiprows = 1)
data0z = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_0_outputs_zv.txt', skiprows = 1)
# Velocities

# RART MED PLOTFIL?

# With interactions

data1x = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_1_outputs_xv.txt', skiprows = 1)
data1y = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_1_outputs_yv.txt', skiprows = 1)
data1z = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_1_outputs_zv.txt', skiprows = 1)

# Plot:
fig, (ax1, ax2, ax3)= plt.subplots(1,3, sharex = True, sharey=True)
fig.suptitle('Phase-space plots w/ and w/o particle interactions') # fjern senere?

ax1.plot(data0x,)






    #axis = ax%d %interactions
    ax1.plot(x1, y1, '-', color=colors[interactions], label=label)
    ax2.plot(x2, y2, '-', color=colors[interactions], label=label)

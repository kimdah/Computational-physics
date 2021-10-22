# Phase space plots - position vs velocity in each direction
import numpy as np
import matplotlib.pyplot as plt

# without interactions
data0x = np.loadtxt('project3/Results/i_100_d_100_p_2_pi_0_outputs_xv.txt')
data0y = np.loadtxt('project3/Results/i_100_d_100_p_2_pi_0_outputs_yv.txt')
data0z = np.loadtxt('project3/Results/i_100_d_100_p_2_pi_0_outputs_zv.txt')

# With interactions
data1x = np.loadtxt('project3/Results/i_100_d_100_p_2_pi_1_outputs_xv.txt')
data1y = np.loadtxt('project3/Results/i_100_d_100_p_2_pi_1_outputs_yv.txt')
data1z = np.loadtxt('project3/Results/i_100_d_100_p_2_pi_1_outputs_zv.txt')

# Plot:
fig, (ax1, ax2)= plt.subplots(1,3, sharex = True, sharey=True)
fig.suptitle('Phase-space plots w/ and w/o particle interactions') # fjern senere?

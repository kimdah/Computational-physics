# relative error  agianst time for 5 different step sizes for 1 particle
import numpy as np
import matplotlib.pyplot as plt


for j in range(0,2):
    if (j == 0):
        method = 'EC'
    else:
        method = 'RK4'

    for i in range(1,6):
        iter = 10**i
        data = np.loadtxt('./Results/%s_i_%d_d_100_p_1_pi_1_outputs_txyz_pert_0_rs_0_f_0.0_w_v_0.0.txt'%(method,iter), skiprows=1)

        t = np.array(data[:,0])
        x = np.array(data[:,1])
        y = np.array(data[:,2])
        z = np.array(data[:,3])

        # Analytical
        x_a =

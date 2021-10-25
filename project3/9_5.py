# relative error  agianst time for 5 different step sizes for 1 particle
# also plot the error convergence rate
import numpy as np
import matplotlib.pyplot as plt

#max_error = np.zeros(5)


for j in range(0,2):
    if (j == 0):
        method = 'EC'
    else:
        method = 'RK4'

    for i in range(1,6):
        iter = 10**i
        data = np.loadtxt('./Results/%s_i_%d_d_100_p_1_pi_1_outputs_txyzv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%(method,iter), skiprows=1)

        t = np.array(data[:,0])
        x = np.array(data[:,1])
        y = np.array(data[:,2])
        z = np.array(data[:,3])
        v = np.array(data[:,4])

        r = np.array(x,y,z)

        # --------- Error convergence rate (9.6)--------
        # max_error[i-1] = np.max(....)




        #------Analytical solution-----
        #Fetches the initial values and propertied for the particle
        x_0 = x[0], y_0 = 0, z_0 = z[0], v_0 = 1    #postion and velocity of particle

        m = 40.078, q=1     #charge and mass

        #Defined constants from task
        B_0 = 9.65*10**1, V_0=9.65*10**8, d=10**4

        #Constants
        omega_0 =(q*B_0)/m , omega_z = np.sqrt((2*q*V_0)/(m*d**2))

        #constants
        omega_plus = (omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2))/2
        omega_minus = (omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2))/2

        #constants
        A_plus = (x_0*omega_minus + v_0)/(omega_minus- omega_plus)
        A_minus = -((v_0+x_0*omega_plus)/(omega_minus- omega_plus))

        #calculates exact values based on t(at timestep) and initial conditions
        x_exact = A_plus*np.cos(omega_plus*t) + A_minus*np.cos(omega_minus*t)
        y_exact = A_plus*np.sin(omega_plus*t) + A_minus*np.sin(omega_minus*t)
        z_exact = z_0*np.cos(omega_z*t)


# x_exact = np.zeros(len(t))
# y_exact = np.zeros(len(t))
# z_exact = np.zeros(len(t))
#
# for i in range(len(t)):
#     x_exact[i] = A_plus*np.cos(omega_plus*t[i]) + A_minus*np.cos(omega_minus*t[i])
#     y_exact[i] = A_plus*np.sin(omega_plus*t[i]) + A_minus*np.sin(omega_minus*t[i])
#     z_exact[i] = z_0*np.cos(omega_z*t[i])

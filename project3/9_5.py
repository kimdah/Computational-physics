# relative error  agianst time for 5 different step sizes for 1 particle
import numpy as np
import matplotlib.pyplot as plt


for j in range(0,2):
    if (j == 0):
        method = 'EC'
    else:
        method = 'RK4'

    for i in range(5,6):
        iterations = 10**i
        data = np.loadtxt('./Results/%s_i_%d_d_100_p_1_pi_1_outputs_txyzv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%(method,iterations), skiprows=1)

        t = np.array(data[:,0])
        x = np.array(data[:,1])
        y = np.array(data[:,2])
        z = np.array(data[:,3])
        v_x = np.array(data[:,4])
        v_y = np.array(data[:,5])
        v_z = np.array(data[:,6])


        #------Analytical solution-----
        #Fetches the initial values and propertied for the particle
        x_0 = x[0]
        y_0 = 0
        z_0 = z[0]
        v_0 = v_y[0]    

        m = 40.078
        q=1     #charge and mass

        #Defined constants from task
        B_0 = 9.65*10**1
        V_0=9.65 *10**8
        d=10**4

        #Constants
        omega_0 =(q*B_0)/m 
        omega_z = np.sqrt((2*q*V_0)/(m*d**2))

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

        relative_error = np.sqrt((x-x_exact)**2+(y-y_exact)**2+(z-z_exact)**2)/np.sqrt((x_exact)**2+(y_exact)**2+(z_exact)**2)
        #print("Max() = ", np.max(np.sqrt((x-x_exact)**2+(y-y_exact)**2+(z-z_exact)**2)))
        #print("Min() = ", np.min(np.sqrt((x_exact)**2+(y_exact)**2+(z_exact)**2)))

        plt.plot(x,y,label="est")
        plt.plot(x_exact,y_exact,label="ex")

        #plt.yscale("log")
        #plt.plot(t,relative_error, label ="iterations= "+str(iterations))

    plt.xlabel("Time in $\mu s$")
    plt.ylabel("Relative error")
    plt.title(method)
    plt.legend()
    plt.grid()

    plt.show()


# x_exact = np.zeros(len(t))
# y_exact = np.zeros(len(t))
# z_exact = np.zeros(len(t))
#
# for i in range(len(t)):
#     x_exact[i] = A_plus*np.cos(omega_plus*t[i]) + A_minus*np.cos(omega_minus*t[i])
#     y_exact[i] = A_plus*np.sin(omega_plus*t[i]) + A_minus*np.sin(omega_minus*t[i])
#     z_exact[i] = z_0*np.cos(omega_z*t[i])

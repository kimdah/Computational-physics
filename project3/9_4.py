from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

#fig = plt.figure()


SMALL_SIZE = 12
MEDIUM_SIZE = 12
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


for i in range(0,2):
    #iter = 10**i
    data = np.loadtxt('./Results/RK4_i_10000_d_100_p_2_pi_%d_outputs_xyz_pert_0_rs_0_f_0.0_w_v_0.0.txt' %i, skiprows=1) #%iter
    x1 = np.array(data[:,0])
    y1 = np.array(data[:,1])
    z1 = np.array(data[:,2])
    x2 = np.array(data[:,3])
    y2 = np.array(data[:,4])
    z2 = np.array(data[:,5])

    # syntax for 3-D projection
    ax = plt.axes(projection ='3d')


    # plotting
    ax.plot3D(x1, y1, z1, label='p1')
    ax.plot3D(x2, y2, z2, label='p2')
    ax.set(xlabel='x(t)', ylabel='y(t)', zlabel='z(t)')
    ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))


    if i==0:
        # no interactions
        inter = ' w/0 interactions'
    else:
        inter = ' w/ interactions'
        # ax.set_xlim([-10**4, 10**4])
        # ax.set_ylim([-10**4, 10**4])
        # ax.set_zlim([-10**4, 10**4])


    ax.set_title('3D plot of the trajectory' + inter)
    #plt.legend()
    plt.savefig('./Figures/xyz%d.pdf'%i)
    plt.show()
    plt.close()

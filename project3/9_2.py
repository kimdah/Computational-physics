# Plot of 2 particles and their motion in xy-plane with and without particle interactions
import numpy as np
import matplotlib.pyplot as plt

# pi:   0: no interactions, 1:interactions

data0 = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_0_outputs_txy.txt', skiprows=1)
data1 = np.loadtxt('./Results/RK4_i_100_d_100_p_2_pi_1_outputs_txy.txt', skiprows=1)

# Particle 1
x01 = np.array(data0[:,1])
y01 = np.array(data0[:,2])
x11 = np.array(data1[:,1])
y11 = np.array(data1[:,2])
# Particle 2
x02 = np.array(data0[:,3])
y02 = np.array(data0[:,4])
x12 = np.array(data1[:,3])
y12 = np.array(data1[:,4])

fig, (ax1, ax2)= plt.subplots(1,2, sharex = True, sharey=True)
fig.suptitle('Motion in xy-plane w/ and w/o particle interactions') # remove later?

# Plotting no interaction vs interaction:

ax1.plot(x01,y01,label='Particle 1')
ax1.plot(x02,y02,label='Particle 2')
ax1.grid()


ax2.plot(x11,y11,label='Particle 1')
ax2.plot(x12,y12,label='Particle 2')
ax2.grid()

ax1.set(xlabel = 'x', ylabel='y', title='No interaction')
ax2.set(xlabel = 'x', title = 'Interaction')

#plt.legend()
plt.legend(bbox_to_anchor=(1.1, 1.05))
plt.savefig('9_1.pdf')
plt.show()

"""
# plotting Particle 1 vs 2:
fig, (ax1, ax2)= plt.subplots(1,2, sharex = True, sharey=True)
fig.suptitle('Motion in xy-plane w/ and w/o particle interactions')

ax1.plot(x01,y01,label='w/o particle interactions')
ax1.plot(x11,y11,label='w/ interactions')

ax2.plot(x02,y02,label='w/o interactions')
ax2.plot(x12,y12,label='w/ interactions')
ax1.set(xlabel = 'x', ylabel='y', title='Particle 1')
ax2.set(xlabel = 'x', title = 'Particle 2')

ax1.grid()
ax2.grid()
plt.legend()
plt.savefig('9_1.pdf')
plt.show()
"""

# Plot of 2 particles and their motion in xy-plane with and without particle interactions
import numpy as np
import matplotlib as plt



fig, (ax1, ax2)= plt.subplots(2, sharey=True) #(ax1, ax2)
fig.suptitle('Motion in xy-plane w/ and w/o particle interactions')



ax1.plot(x1, y1, '-', label='Particle1 w/inter.')
ax1.plot(x1o, y1o, '-', color='r', label='Particle1 w/o inter.')

ax1.plot(x2, y2, '-', label='Particle2 w/inter.')
ax1.plot(x2o, y2o, '-', color='r', label='Particle2 w/o inter.')


ax1.set_title('Particle 1')
ax2.set_title('Particle 2')

ax1.set(xlabel = 'x')
ax2.set(xlabel = 'x')

# LEGENDS!

plt.ylabel('y')
plt.savefig('9_2plot.pdf')
plt.grid()
plt.show()

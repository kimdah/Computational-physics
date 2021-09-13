import matplotlib.pyplot as plt
import numpy as np
import math

k = 3 # max exponent of 10 number of points
colors = ['r', 'g', 'y'] # add more

fig, (ax1, ax2) = plt.subplots(2, sharex=True)
fig.suptitle('Errors')

for i in range(1, k+1):
    N = 10**i
    errors = np.loadtxt('errors%d.txt' %N) # only once?
    x = np.array(errors[:,0])
    abs_error = np.array(errors[:,1])
    rel_error = np.array(errors[:,2])

    #plt.plot(x, abs_error, 'o-', color=colors[i-1], label='Approximation v*(x) with n =%d' %n)
    ax1.plot(x, abs_error, '-', color=colors[i-1], label='Abs_error with %d' %N)
    ax2.plot(x, rel_error, '-', color=colors[i-1], label='Rel_error with %d' %N)




#plt.xlabel('x')
#plt.ylabel('u(x)') # v? !!!!
#plt.title('General algorithm approximation vs exact of Poisson eq.')
#plt.plot(x, u, color='k', label='Exact u(x)') # Plotting u(x) with highest n
plt.legend()
plt.savefig('error_plot.pdf')
plt.show()











"""
plt.xlabel('x')
plt.ylabel('log10(|u(x)-v(x)|)')

plt.plot(x, abs_error, color='r', label='n = ' + npoints)
plt.plot(x, abs_error, color='g', label='n = ' + npoints)



plt.savefig('abs_error.pdf')


for ax in axs.flat:
    ax.set(xlabel='x-label', ylabel='y-label')

"""

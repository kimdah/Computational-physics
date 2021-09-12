import matplotlib.pyplot as plt
import numpy as np
import math

k = 3 # max exponent of 10 number of points
colors = ['r', 'g', 'y'] # add more

for i in range(1, k+1):
    n = 10**i
    exact_data = np.loadtxt('exact_data%d.txt' %n) # only once?
    x = np.array(exact_data[:,0])
    u = np.array(exact_data[:,1])
    genapprox_data = np.loadtxt('approx_general%d.txt' %n)
    v = np.array(genapprox_data[:,1])

    vstar = np.array([0]) # boundary point u_0 = 0
    vstar = np.append(vstar, v)
    vstar = np.append(vstar, 0) # appending boundary point u_1 = 0


    
    for (i : abs_error):
        abs_error = math.log10(abs(u-vstar))
        rel_error = math.log10(abs(u-vstar)/u)

    plt.plot(x, abs_error, 'o-', color=colors[i-1], label='Approximation v*(x) with n =%d' %n)


plt.xlabel('x')
plt.ylabel('u(x)') # v? !!!!
plt.title('General algorithm approximation vs exact of Poisson eq.')
plt.plot(x, u, color='k', label='Exact u(x)') # Plotting u(x) with highest n
plt.legend()
plt.savefig('approx_general.pdf')
plt.show()












plt.xlabel('x')
plt.ylabel('log10(|u(x)-v(x)|)')

plt.plot(x, abs_error, color='r', label='n = ' + npoints)
plt.plot(x, abs_error, color='g', label='n = ' + npoints)



plt.savefig('abs_error.pdf')


for ax in axs.flat:
    ax.set(xlabel='x-label', ylabel='y-label')

import matplotlib.pyplot as plt
import numpy as np
import math

# Plotting for multiple n values
k = 3 # max exponent of 10 number of points
colors = ['b', 'r', 'g', 'y'] # add more

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
    plt.plot(x, vstar, 'o', color=colors[i], label='Approximation v*(x) with n =%d' %n)


plt.xlabel('x')
plt.ylabel('u(x)') # v? !!!!
plt.title('General algorithm approximation vs exact of Poisson eq.')
plt.plot(x, u, color='k', label='Exact u(x)') # Plotting u(x) with highest n
plt.legend()
plt.savefig('approx_general.pdf')

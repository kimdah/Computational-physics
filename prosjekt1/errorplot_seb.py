import matplotlib.pyplot as plt
import numpy as np
import math

k=3
colors = ['b', 'r', 'g', 'y']

for i in range (1, k+1):
    n = 10**i
    exact_data = np.loadtxt('exact_data%d.txt' %n) # only once?
    x = np.array(exact_data[:,0])
    u = np.array(exact_data[:,1])
    genapprox_data = np.loadtxt('approx_general%d.txt' %n)
    v = np.array(genapprox_data[:,1])

    vstar = np.array([0]) # boundary point u_0 = 0
    vstar = np.append(vstar, v)
    vstar = np.append(vstar, 0) # appending boundary point u_1 = 0
    np.log10(abs_error) = np.log10(np.abs(np.array(u)-np.array(vstar)))
    plt.plot(x, abs_error, color=colors[i], label='n = %d' %n)

plt.xlabel('x')
plt.ylabel('log10(|u(x)-v(x)|)')


#plt.plot(x, abs_error, color='g', label='n = ' + n)


#plt.show()
plt.savefig('abs_error.pdf')

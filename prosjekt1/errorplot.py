import matplotlib.pyplot as plt
import numpy as np
import math

data = np.loadtxt('data.txt')
x = data[:,0]
u = data[:,1]
datan = np.loadtxt('datan.txt')
v = datan[:,0] //??

abs_error = math.log10(abs(u-v))




plt.xlabel('x')
plt.ylabel('log10(|u(x)-v(x)|)')

plt.plot(x, abs_error, color='r', label='n = ' + npoints)
plt.plot(x, abs_error, color='g', label='n = ' + npoints)



plt.savefig('abs_error.pdf')


for ax in axs.flat:
    ax.set(xlabel='x-label', ylabel='y-label')

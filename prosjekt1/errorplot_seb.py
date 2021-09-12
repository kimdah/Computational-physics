import matplotlib.pyplot as plt
import numpy as np
import math

exact_data = np.loadtxt('exact_data.txt') # only once?
x = np.array(exact_data[:,0])
u = np.array(exact_data[:,1])
genapprox_data = np.loadtxt('approx_general.txt')
v = np.array(genapprox_data[:,1])

vstar = np.array([0]) # boundary point u_0 = 0
vstar = np.append(vstar, v)
vstar = np.append(vstar, 0) # appending boundary point u_1 = 0

abs_error = math.log10(np.absolute(np.array(u)-np.array(vstar)))

plt.xlabel('x')
plt.ylabel('log10(|u(x)-v(x)|)')

plt.plot(x, abs_error, color='r', label='n = ' + n)
#plt.plot(x, abs_error, color='g', label='n = ' + n)


plt.show()
#plt.savefig('abs_error.pdf')

import matplotlib.pyplot as plt
import numpy as np
import math

k=3
colors = ['b', 'r', 'g', 'y']

for i in range (1, k+1):
    n = 10**i
    exact_data = np.loadtxt('exact_data%d.txt' %n)
    x = np.array(exact_data[:,0])
    u = np.array(exact_data[:,1])
    genapprox_data = np.loadtxt('approx_general%d.txt' %n)
    v = np.array(genapprox_data[:,1])
    # Prover aa ekskludere boundary:
    print(x.size, v.size)
    #print(x)
    x = np.delete(x, [0, x.size-1])
    print(x)
    #print(u)
    u = np.delete(u, [0, u.size-1])
    #print(u)
    print(x.size, u.size, v.size)

    abs_err = np.abs(u - v) # Dangerous for loss of precision!!!!!!!!
    rel_err = np.abs(abs_err / u)

    # Ask it to ignore results that are zero due to v-u = 0, and log 0 being undef.
    log10_abs_err = np.log10(abs_err)#, out=abs_err, where=abs_err > 0)
    log10_rel_err = np.log10(rel_err)#, out=rel_err, where=rel_err > 0)
    plt.plot(x, log10_abs_err, color=colors[i], label='n = %d' %n)

plt.title('Absolute error')
plt.xlabel('x')
plt.ylabel('log10(|u(x)-v(x)|)')


#plt.plot(x, abs_error, color='g', label='n = ' + n)


#plt.show()
plt.savefig('abs_error.pdf')

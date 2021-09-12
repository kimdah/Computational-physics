
#from script import *
import matplotlib.pyplot as plt
import numpy as np
data = np.loadtxt('exact_data.txt')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.plot(data[:,0], data [:,1], 'ro')
plt.savefig('problem2.pdf')

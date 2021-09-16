import matplotlib.pyplot as plt
import numpy as np
data = np.loadtxt('exact_data1000.txt')
plt.title('Exact solution of Poisson eq.\nHere with N=1000')
plt.xlabel('x')
plt.ylabel('u(x)= 1- (1-exp(-10)x - exp(-10*x))')
plt.plot(data[:,0], data [:,1])
plt.savefig('problem2.pdf')

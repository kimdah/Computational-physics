import matplotlib.pyplot as plt
import numpy as np

colors = ['r', 'g', 'b', 'r', 'c', 'm', 'y'] # shorten later

for n in range (10, 100, )(n=10; n<= 100; n*10){


}
for i in range(1, 3):
    n = 10**i # n = 10, 100
    data = np.loadtxt('output%d.txt' %n)
    xscaled = np.array(data[:,0])
    eigvec1 = np.array(data[:,1])
    eigvec2 = np.array(data[: 2])
    eigvec3 = np.array(data[:,3])
    eigval = np.array(np.array(data[:,1]),np.array(data[: 2]), np.array(data[:,3])) # isteden?

    xhat = np.array();
    xhat.append(0); # boundary value xhat = 0
    xhat = np.append(xhat, xscaled)
    xhat.append(1) # boundary value xhat = 1

    for j in range(1,4):

        plt.plot(xhat, vstar, color=colors[i-1], label='n =%d' %n)



    plt.ylabel('Eigenvectors')
    plt.xlabel('xhat')
    plt.savefig('x_vs_eigvecn=%d.pdf' %n)
    plt.show()

import matplotlib.pyplot as plt
from math import *
import numpy as np
import sys



#defines x, u and v lists
N =[]
itterations =[]


#----------------Task2: plotting u------------
with open("task6_dataset.txt") as myfile: #opens file
    myfile.readline()
    for line in myfile:                 #reads each line of the file
        mydata = line.split(",")       #splits each line with the seperator
        N.append(float(mydata[0]))      #Appends values to x list
        itterations.append(float(mydata[1]))      #Appends values to u list

myfile.close()  #closes file


N = np.array(N)
itterations = np.array(itterations)


#plots function u(x) 

plt.plot(N,itterations, label ="Number of rotations")   




plt.xlabel("N values")
plt.ylabel("Number of rotations")
plt.title("Number of rotations when A is a symmetrical tridiagonal NxN matrix")
plt.legend()
plt.grid()     
plt.savefig('Task6plot.pdf')
plt.show()      


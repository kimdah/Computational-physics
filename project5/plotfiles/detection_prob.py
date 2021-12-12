import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
import sys

#Adjusting text size
SMALL_SIZE = 14
MEDIUM_SIZE = 14
# BIGGER_SIZE = 17
#
# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

#Input arguments that determines what is being plotted.
filename = sys.argv[1]
slits = sys.argv[2]
#x_axis_label = sys.argv[2] # if we want general one
#y_axis_label = sys.argv[3]

slice = pa.cx_mat() #Create pa.mat object (just as arma::mat in C++)
slice.load("./datafiles/"+str(filename)) #Load the content of the matrix you saved into your Python program.
slice = np.array(slice)

h = 0.005
points = int((1/h) +1) # M = 201
y = np.linspace(0,1,points)

index = int(0.8/h) # x = 0.8 / stepsize h
U = slice[index,:]
p = np.real(np.conj(U)*U); # here???? Shouldnt this be done to the entire matrix?
p_norm = p/np.sqrt(np.sum(p**2))# normalise
#print("Sum prob after normalising: " , np.sum(p_norm))

plt.plot(y, p_norm)
plt.xlabel("y [Units of distance/1]")
plt.ylabel("p(y|x=0.8; t=0.002)")
plt.grid()

plt.subplots_adjust(
top=0.915,
bottom=0.165,
left=0.20,
right=0.95,
hspace=0.0,
wspace=0.0
)

plt.savefig("./figures/detection_prob_slits_%s.pdf" %slits)
#plt.show()

# Studying the burn-in time for L = 20
import numpy as np
from plotlib import makeplots
import matplotlib.pyplot as plt







#Adjusting text size
SMALL_SIZE = 17
MEDIUM_SIZE = 17
BIGGER_SIZE = 17

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


plotfiles = ['test.txt','test2.txt']	#Name of file
x_name = "$\lambda$"					#xlabel
y_name = "$\eta$"						#ylabel

#Use 'seperate' to save as different plots and 'same' to plot in same plot

makeplots(plotfiles, x_name, y_name, 'seperate')
makeplots(plotfiles, x_name, y_name, 'same')

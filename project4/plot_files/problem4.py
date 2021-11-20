# Studying the burn-in time for L = 20
import numpy as np
#from plotlib import makeplots
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


data_values =  ['samples','E','M', 'expval_eps', 'expval_m', 'C_V', 'sucept']
data = np.loadtxt('../datafiles/task4.txt',skiprows=1)
samples = np.array(data[:,0])
E = np.array(data[:,1])
M = np.array(data[:,2])
expval_eps = np.array(data[:,3])
expval_m = np.array(data[:,4])
C_V = np.array(data[:,5])
sucept = np.array(data[:,6])



fig, ax = plt.subplots(figsize = (6, 5))
plt.subplots_adjust(
top=0.95,
bottom=0.15,
left=0.15,
right=0.98,
hspace=0.2,
wspace=0.2
)

for i in range(3,7):

	plt.plot(np.array(data[:,0]),np.array(data[:,i]), label=data_values[i])
	#ax.set_title("Training mse")
	ax.set_ylabel(data_values[i])
	ax.set_xlabel(data_values[0])
	plt.grid()
	plot_name = "Task4_"+data_values[0]+"_"+data_values[i]+"_figure.pdf"
	plt.legend()
	plt.savefig('../figures/'+plot_name)
	plt.clf()

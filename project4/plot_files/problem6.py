# Estimating the pdf of epsilon
import numpy as np
#from plotlib import makeplots
import matplotlib.pyplot as plt
#%matplotlib inline

#Adjusting text size
SMALL_SIZE = 14
MEDIUM_SIZE = 17
BIGGER_SIZE = 17

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

fig, ax = plt.subplots(figsize = (6, 5))
plt.subplots_adjust(
top=0.95,
bottom=0.15,
left=0.2,
right=0.95,
hspace=0.2,
wspace=0.2
)
#ax.set_ylabel(data_units[i])
#ax.set_xlabel(data_units[0])
#plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

# np.random.seed(42)
# x = np.random.normal(size=1000)


# CHOSE UNORDERED - MAYBE DO BOTH?
# dataT1u = np.loadtxt('./datafiles/ncyc_1e5_L_20_T_1.0_unordered.txt', skiprows=1)
# dataT2u = np.loadtxt('./datafiles/ncyc_1e5_L_20_T_2.4_unordered.txt', skiprows=1)
#
# eps1 = np.array(dataT1u[:,3])
# eps2 = np.array(dataT2u[:,3])

for i in range(0,1): # 2
    if i == 0:
        T = 1.0
    else:
        T = 2.4

    data = np.loadtxt('./datafiles/ncyc_1e4_L_20_T_%.1f_unordered.txt' %T, skiprows=1)
    eps= np.array(data[1:,3]) # np.array?
    print(eps[0])
    #eps = eps_temp[]
    print(len(eps))
    # Using Freedman–Diaconis rule to be more scientific in choosing the "right" bin width
    q25, q75 = np.nanpercentile(eps, [0.25, 0.75])
    print(q25, q75)
    bin_width = 2 * (abs(q75) - abs(q25)) * len(eps) ** (-1/3)
    print(bin_width)
    print(max(eps), eps.min())
    bins = round((eps.max() - eps.min()) / bin_width)
    print("Freedman–Diaconis number of bins:", bins)
    #bins = 30

    plt.style.use('seaborn-white')
    #plt.hist(eps, density=True, bins=bins, stacked=True) #density=True?

    #plt.hist(x, density=True, bins=30)  # density=False would make counts
    plt.ylabel('Probability')
    plt.xlabel('Data')
    plt.show()

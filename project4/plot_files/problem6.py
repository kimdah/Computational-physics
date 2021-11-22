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

# CHOSE UNORDERED - MAYBE DO BOTH?
for i in range(0,2): 
    if i == 0:
        T = 1.0
    else:
        T = 2.4

    data = np.loadtxt('./datafiles/ncyc_1e4_L_20_T_%.1f_unordered.txt' %T, skiprows=1)
    eps = np.array(data[1:,3]) # np.array?
    # Using Freedman–Diaconis rule to be more scientific in choosing the "right" bin width
    q25, q75 = np.percentile(eps, [0.25, 0.75])
    #print(q25, q75)
    bin_width = 2 * (q75 - q25) * len(eps) ** (-1/3)
    # print(bin_width)
    # print(max(eps), max(abs(eps)), abs(max(eps)))
    # print(min(eps), min(abs(eps)), abs(min(eps)))
    # print(max(abs(eps))- min(abs(eps)))
    #bin_width = 0.01
    bins = round((max(abs(eps))- min(abs(eps))) / bin_width)
    print("Freedman–Diaconis number of bins:", bins)
    #bins = 30
    bins = round(np.sqrt(len(eps))) # 100
    print(np.sqrt(len(eps)), len(eps))
    bins = 32
    print("bin_width: ", round((max(abs(eps))- min(abs(eps)))/bins))

    #plt.style.use('seaborn-white')
    plt.hist(eps, density=True, bins=bins, stacked=True)

    plt.ylabel('$p_{\epsilon}(\epsilon)_{est}$')
    plt.xlabel('$\epsilon$')
    plt.grid()
    plt.savefig('figures/probability_distribution_T='+str(T)+'.pdf')
    plt.clf()
    plt.cla()

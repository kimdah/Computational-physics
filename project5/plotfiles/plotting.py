import pyarma as pa
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

sns.set()

A = pa.cx_mat() #Create pa.mat object (just as arma::mat in C++)
A.load("./datafiles/test.dat") #Load the content of the matrix you saved into your Python program.

fig, ax = plt.subplots(figsize=(5, 5))
sns.heatmap(np.real(A), annot=True, ax=ax, cmap="viridis")

print(A)

plt.show()
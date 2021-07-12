# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data for setup
M = 1.0
mu = 0.5
r = 30
symmetry = 2
# make the plot
fig = plt.figure()
N = 0

# volume integral dataset out
data1 = np.loadtxt("c4ism1/ProcaDensities.dat")
timedata = data1[:,0]
dM = symmetry*data1[:,1] - symmetry*data1[N,1]
plt.plot(timedata, dM, '-', lw = 1.0, label="M-M0 LR")

data1 = np.loadtxt("c4ism1HR/ProcaDensities.dat")
timedata = data1[:,0]
dM = symmetry*data1[:,1] - symmetry*data1[N,1]
plt.plot(timedata, dM, '-', lw = 1.0, label="M-M0 MR")

data1 = np.loadtxt("c4ism1HHR/ProcaDensities.dat")
timedata = data1[:,0]
dM = symmetry*data1[:,1] - symmetry*data1[N,1]
plt.plot(timedata, dM, '-', lw = 1.0, label="M-M0 HR")

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Change in Cloud Mass")
#plt.xlim(0, 1000)
#plt.ylim(1e-1, 1e2)
plt.legend(loc=0)
plt.grid()

# save as png image
filename = "MvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)

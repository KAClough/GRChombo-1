# A simple python script to plot the GW
# signals over time, for a chosen mode

import numpy as np;
import matplotlib.pyplot as plt;

# output data for setup
M = 1.0
mu = 0.5
r = 30
a = 0.99
symmetry = 2
N = 4000 # when to normalise to
alpha = M * mu
r_plus = M + np.sqrt(M*M - a*a)
omega_Re = mu * (1.0 - 0.5 * alpha**2 - 0.33 * alpha **3.0)
omega_Im = 0.99 * alpha**7.0 * (a - 2 * r_plus * omega_Re)

# make the plot
fig = plt.figure()

# volume integral dataset out
data1 = np.loadtxt("c4is0/ProcaDensities.dat")
labelstring = "M/M0 LR"
timedata = data1[:,0]
dM = symmetry*data1[:,1]/(symmetry*data1[int(N/2),1])
plt.semilogy(timedata, dM, '-', lw = 1.0, label=labelstring)

# volume integral dataset out
data1 = np.loadtxt("c4is0HR/ProcaDensities.dat")
labelstring = "M/M0 MR"
timedata = data1[:,0]
dMMR = np.abs(symmetry*data1[:,1]/(symmetry*data1[int(N),1]))
plt.semilogy(timedata, dMMR, '-', lw = 1.0, label=labelstring)

# analytic
plt.semilogy(timedata, dMMR[N]*np.exp(2*omega_Im * (timedata - timedata[N])), '--', lw = 1.0, label="analytic")
plt.semilogy(timedata, dMMR[N]*np.exp(2*0.95*omega_Im * (timedata - timedata[N])), '--', lw = 1.0, label="analytic +")
plt.semilogy(timedata, dMMR[N]*np.exp(2*1.05*omega_Im * (timedata - timedata[N])), '--', lw = 1.0, label="analytic -")

# volume integral dataset out
data1 = np.loadtxt("c4is0HHR/ProcaDensities.dat")
labelstring = "M/M0 HR"
timedata = data1[:,0]
dMHR = np.abs(symmetry*data1[:,1]/(symmetry*data1[N,1]))
plt.semilogy(timedata, dMHR, '-', lw = 1.0, label=labelstring)

# make the plot look nice
plt.xlabel("time")
plt.ylabel("Cloud E/J")
plt.xlim(0, 4000)
plt.ylim(1e-1, 1e1)
plt.legend()
plt.grid()

# save as png image
filename = "EJvsT" + "_mu" + str(mu) + ".png"
plt.savefig(filename)

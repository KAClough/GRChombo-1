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
data1 = np.loadtxt("ProcaDensities.dat")
timedata = data1[:,0]
dM = symmetry*data1[:,1] - symmetry*data1[N,1]
#print(dM)

# flux dataset out
data1 = np.loadtxt("Flux_integrals.dat")
labelstring = "integral(Flux * dt)"
timedata = data1[:,0]
dt = timedata[1] - timedata[0]
NetFlux = data1[:,1] - data1[:,3]
FluxOuter = data1[:,3]
FluxInner = data1[:,1]
Fdt = np.zeros_like(timedata)
FIdt = np.zeros_like(timedata)
FOdt = np.zeros_like(timedata)
for i, F in enumerate(NetFlux) :
    if (i > N) :
        Fdt[i] += Fdt[i-1] + F * dt
        FOdt[i] += FOdt[i-1] + FluxOuter[i] * dt
        FIdt[i] += FIdt[i-1] + FluxInner[i] * dt
    else :
        dM[i] = 0.0

plt.plot(timedata, Fdt, '-', lw = 1.0, label="net flux")
plt.plot(timedata, FIdt, '-', lw = 1.0, label="inner flux")
plt.plot(timedata, FOdt, '-', lw = 1.0, label="outer flux")
labelstring = "M-M0"
plt.plot(timedata, dM, '-', lw = 1.0, label=labelstring)
#print(Fdt)

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

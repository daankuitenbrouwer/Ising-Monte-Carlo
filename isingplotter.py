import numpy as np
import matplotlib.pyplot as plt

# runs and Trange have to be specified to open the mentioned datafile

runs = 100
Trange = np.linspace(1,3,300)
magarray = np.loadtxt("magnetizationruns%sTrange%s.txt"%(runs,len(Trange)))

# Critical temperature and exponent are defined, the theoretical plot and empirical plot are both made.
Tc = 2.2
beta = 0.5
mag = (abs(Trange - Tc))**beta
plt.plot(Trange,mag)
plt.scatter(Trange,np.mean(magarray,axis=0))
plt.show()

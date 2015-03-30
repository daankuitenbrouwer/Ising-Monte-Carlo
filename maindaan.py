import numpy as np
import functionsdaan,isingfunctionsdaan
import matplotlib.pyplot as plt

size = 22 #Width of the lattice. 
init_M = 1 #initial magnetization
Trange = np.linspace(0.1,3,30)
runs = 20
magarray = np.zeros((runs,len(Trange)))

for run in range(runs):
  lattice,BlackList = isingfunctionsdaan.InitializeLattice(size,init_M)
  Magnetization_func_of_T = isingfunctionsdaan.CritTemp(lattice,BlackList,Trange)
  magarray[run] = Magnetization_func_of_T
  print 100*(float(run)/runs),'%'
  
np.savetxt("magnetizationruns%sTrange%s.txt"%(runs,len(Trange)),magarray)
plt.scatter(Trange,np.mean(magarray,axis=0))
plt.show()

import numpy as np
import functionsdaan,isingfunctionsdaan
import matplotlib.pyplot as plt

size = 22 #Width of the lattice. 
init_M = 1 #initial magnetization
lowT = 1
highT = 3
Trange = np.linspace(lowT,highT,400)
runs = 100
magarray = np.zeros((runs,len(Trange)))
maxrecursion = 99999999999

for run in range(runs):
  lattice,BlackList = isingfunctionsdaan.InitializeLattice(size,init_M)
  Magnetization_func_of_T = isingfunctionsdaan.CritTemp(lattice,BlackList,Trange,maxrecursion)
  magarray[run] = Magnetization_func_of_T
  print 100*(float(run)/runs),'%'
  
np.savetxt('magnetizationruns%sTsteps%ssize%slowT%shighT%s.txt'%(runs,len(Trange),size,lowT,highT),magarray)
plt.scatter(Trange,np.mean(magarray,axis=0))
plt.show()


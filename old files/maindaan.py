import numpy as np
import functionsdaan,isingfunctionsdaan
import matplotlib.pyplot as plt
import sys

sys.setrecursionlimit(100000)

size = 100 #Width of the lattice. 
init_M = 1 #initial magnetization
lowT = 0                               
highT = 2.5
Trange = np.linspace(lowT,highT,100)
runs = 1
magarray = np.zeros((runs,len(Trange)))

for run in range(runs):
  lattice,BlackList = isingfunctionsdaan.InitializeLattice(size,init_M)
  Magnetization_func_of_T = isingfunctionsdaan.CritTemp(lattice,BlackList,Trange)
  magarray[run] = Magnetization_func_of_T
  print 100*(float(run)/runs),'%'
  
np.savetxt('magnetizationruns%sTsteps%ssize%slowT%shighT%s.txt'%(runs,len(Trange),size,lowT,highT),magarray)
plt.scatter(Trange,np.mean(abs(magarray),axis=0))
plt.show()


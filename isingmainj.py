import numpy as np
import isingfunctionsj
import matplotlib.pyplot as plt
import sys

sys.setrecursionlimit(100000)

# settings
init_M = 1 #initial magnetization
runs = 10
Tsteps = 200
size = 100 #Width of the lattice. 
lowT = 1.3                              
highT = 4.7
Tc = 2.2727
Trange = np.linspace(lowT,highT,Tsteps)
lowh = 0
highh = 0.2
hsteps = 100
hrange = np.linspace(lowh,highh,hsteps)
Niter = 1000 # number of new random position to let the system find an equilibrium

magarray = np.zeros((runs,len(Trange)))
en_array = np.zeros((runs,len(Trange)))
mag_h_array = np.zeros((runs,len(hrange)))


'''
for run in range(runs):
  lattice,BlackList = isingfunctionsj.InitializeLattice(size,init_M)
  Magnetization_func_of_T,energy_func_of_T = isingfunctionsj.CritTemp(lattice,BlackList,Trange,Niter,lowh)
  magarray[run] = Magnetization_func_of_T
  en_array[run] = energy_func_of_T
  print 100*(float(run)/runs),'%'
  
  


  
np.savetxt('magnetizationruns%sTsteps%ssize%slowT%shighT%s.txt'%(runs,len(Trange),size,lowT,highT),magarray)
np.savetxt('energyruns%sTsteps%ssize%slowT%shighT%s.txt'%(runs,len(Trange),size,lowT,highT),en_array)
plt.scatter(Trange,np.mean(abs(magarray),axis=0))
#plt.scatter(Trange,np.mean(en_array,axis=0))
#plt.scatter(Trange,(-0.5*4*len(lattice)**2)*np.mean(abs(magarray),axis=0)**2)
#

plt.show()
#'''

for run in range(runs):
  lattice,BlackList = isingfunctionsj.InitializeLattice(size,0)
  Mag_func_of_h = isingfunctionsj.Crit_Exp_Delta(lattice,BlackList,hrange,Niter,Tc)
  mag_h_array[run] = Mag_func_of_h
  print 100*(float(run)/runs),'%'
  
np.savetxt('magnetizationruns%shsteps%ssize%slowh%shighh%s.txt'%(runs,hsteps,size,lowh,highh),mag_h_array)
  
plt.scatter(hrange,abs(np.mean(mag_h_array,axis=0)))
plt.show()


  
  
  

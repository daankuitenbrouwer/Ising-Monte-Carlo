import matplotlib.pyplot as plt
import numpy as np
import scipy
import random

import initising
import ising_functions as isf


#Initialize, still to be put in a class somewhere else
N = 100
Nstep = 100000
T = 1 #taking T smaller than 0.1 may yield infinite values.
Kb = 1
Beta = 1/(T*Kb)
magfield = 0.1

#Initialize
Lattice,marray,parray,dearray = initising.init_Lattice(N,Nstep)
initialenergy = initising.initial_energy(Lattice,N)

#print Lattice, marray



#the loop that flips spins
cnt = 1
for step in range(Nstep):
  x = np.random.randint(0,N) #positions of particle in lattice that is to be swapped
  y = np.random.randint(0,N)
  #print x,y
  formerlattice = np.empty_like(Lattice) #deep copy
  formerlattice[:] = Lattice 
  Lattice[y][x] = -Lattice[y][x]
  de,dm = isf.calcdedm(Lattice,Nstep,magfield,x,y,N)
  p = np.exp(de*Beta)
  parray[step] = p
  dearray[step] = de
  #print p, de, de*Beta, dm
  take = isf.accdeclnewlat(p)
  #print marray[step - 1]
  if take == 0:
    Lattice = formerlattice
    if step == 0:
      pass
    else:
      marray[step] = marray[step - 1]
  else:
    if step == 0:
      pass
    else:
      marray[step] = marray[step - 1] + dm
    
#print marray
    
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(range(len(marray)),(marray/(N**2)))
#ax1.plot(range(len(parray)),parray)
#ax1.plot(range(len(dearray)),dearray)
plt.show()



from __future__ import division
from scipy.ndimage import measurements as sci
import numpy as np
import functions
import functionsdaan
#import ani
import matplotlib.pyplot as plt


size =30 #Width of the lattice. The lattice should always be square, otherwise problems occur in de periodic bc definition.
init_M = 1 #initial magnetization
lattice = functionsdaan.initializelat(size,init_M) #initialize the girl, M sets magnetization
while np.average(lattice) != init_M:
    lattice =functionsdaan.initializelat(size,init_M)
blacklist,futurelattice = functionsdaan.initializeelse(size)
T = 0.8  #temperature, entropy vs potential energy
invT = 1./T
beta = invT
runtime = 5000 #Total time of the simulation
B = 0.1 #The favorite orientation
Trange = np.linspace(0.001,1,50)

#Initializing arrays for clusterflippin'
xval,yval = functionsdaan.rand_init_pos(blacklist,lattice) 
blacklist[yval][xval] = 1 # I flip the first spin unconditionally
lattice[yval][xval] = -1*lattice[yval][xval] #the module about accepting or rejecting the first randomly picked spin must be built still
futurelattice[yval][xval] = lattice[yval][xval]
comp_xval = xval + 1
comp_yval = yval
cnt = 0


magnetization = np.zeros(len(Trange))
for T in Trange:
  beta = 1./T
  futurelattice = functionsdaan.wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)
  magnetization[T] = np.average(futurelattice)
  futurelattice[:] = -1*futurelattice[:]
  print T/len(Trange)
  
plt.plot(Trange,magnetization)
plt.show()


print magnetization

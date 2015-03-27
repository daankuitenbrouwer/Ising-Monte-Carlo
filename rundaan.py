from __future__ import division
from scipy.ndimage import measurements as sci
import numpy as np
import functions
import functionsdaan
#import ani
import matplotlib.pyplot as plt


size =1000 #Width of the lattice. The lattice should always be square, otherwise problems occur in de periodic bc definition.
init_M = 1 #initial magnetization
lattice = functionsdaan.initializelat(size,init_M) #initialize the girl, M sets magnetization
while np.average(lattice) != init_M:
    lattice =functionsdaan.initializelat(size,init_M)
blacklist,futurelattice = functionsdaan.initializeelse(size)
runs = 100
Trange = np.linspace(0.1,3,100)
magtarray=np.zeros((runs,len(Trange)))


#Initializing arrays for clusterflippin'
xval,yval = functionsdaan.rand_init_pos(blacklist,lattice) 
blacklist[yval][xval] = 1 # I flip the first spin unconditionally
lattice[yval][xval] = -1*lattice[yval][xval] #the module about accepting or rejecting the first randomly picked spin must be built still
futurelattice[yval][xval] = lattice[yval][xval]
comp_xval = xval + 1
comp_yval = yval


for ind,run in enumerate(range(runs)):
  magtarray[run][:] = functionsdaan.crittemp(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,Trange,size)
  print ind/runs,'%'
print magtarray 
np.savetxt("magnetizationruns%sTrange%s.txt"%(runs,len(Trange)),magtarray)
plt.scatter(Trange,np.mean(magtarray,axis=0))
plt.show()
  



from __future__ import division
from scipy.ndimage import measurements as sci
import numpy as np
import functions
import functionsdaan
#import ani
import matplotlib.pyplot as plt


size =5 #Width of the lattice. The lattice should always be square, otherwise problems occur in de periodic bc definition.
init_M = 1 #initial magnetization
lattice = functionsdaan.initializelat(size,init_M) #initialize the girl, M sets magnetization
while np.average(lattice) != init_M:
    lattice =functionsdaan.initializelat(size,init_M)
blacklist,futurelattice = functionsdaan.initializeelse(size)
#T = 0.3  #temperature, entropy vs potential energy
#invT = 1./T
#beta = invT
runtime = 5000 #Total time of the simulation
B = 0.1 #The favorite orientation
Trange = np.linspace(0.1,1,10)

#Initializing arrays for clusterflippin'
xval,yval = functionsdaan.rand_init_pos(blacklist,lattice) 
blacklist[yval][xval] = 1 # I flip the first spin unconditionally
lattice[yval][xval] = -1*lattice[yval][xval] #the module about accepting or rejecting the first randomly picked spin must be built still
futurelattice[yval][xval] = lattice[yval][xval]
comp_xval = xval + 1
comp_yval = yval
cnt = 0

beta = 1./Trange[0] # test start
#plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
#plt.show()

futurelattice = functionsdaan.wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)
'''
plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
plt.show()

futurelattice = functionsdaan.wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)

plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
plt.show()
 '''
 
print Trange

magnetization = np.zeros(len(Trange))
for ind,T in enumerate(Trange):
  beta = 1./T
  #plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
  #print 'fut'
  #plt.show()
  
  lattice[:] = futurelattice[:]
  #plt.imshow(lattice,cmap="Greys",interpolation="nearest")
  #print 'lat'
  #plt.show()
  blacklist,futurelattice = functionsdaan.initializeelse(size)
  xval,yval = functionsdaan.rand_init_pos(blacklist,lattice) 
  blacklist[yval][xval] = 1 # I flip the first spin unconditionally
  lattice[yval][xval] = -1*lattice[yval][xval] #the module about accepting or rejecting the first randomly picked spin must be built still
  futurelattice[yval][xval] = lattice[yval][xval]
  comp_xval = xval + 1
  comp_yval = yval
  xdir = comp_xval - xval
  ydir = comp_yval - yval
  xdir,ydir,xval,yval,comp_xval,comp_yval = functionsdaan.periodic_boundary(lattice,xdir,ydir,xval,yval,comp_xval,comp_yval)
  cnt = 0
  futurelattice = functionsdaan.wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)
  magnetization[ind] = np.average(futurelattice)
  futurelattice[:] = -1*futurelattice[:]
  print magnetization[ind],'magnetization[ind]'
  #plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
  #plt.show()
plt.plot(Trange,magnetization)
plt.show()


print magnetization


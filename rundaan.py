from __future__ import division
from scipy.ndimage import measurements as sci
import numpy as np
import functions
import functionsdaan
#import ani
import matplotlib.pyplot as plt


size =6 #Width of the lattice. The lattice should always be square, otherwise problems occur in de periodic bc definition.
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
Trange = np.linspace(0.001,1,5)

#Initializing arrays for clusterflippin'
xval = 1 #np.random.randint(0,len(lattice)+1)
yval = 1 #np.random.randint(0,len(lattice)+1)
blacklist[yval][xval] = 1
lattice[yval][xval] = -1*lattice[yval][xval] #the module about accepting or rejecting the first randomly picked spin must be built still
futurelattice[yval][xval] = lattice[yval][xval]
comp_xval = xval + 1
comp_yval = yval
blacklist[comp_yval][comp_xval] = 1
cnt = 0
print len(lattice)
futurelattice = functionsdaan.wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)
print futurelattice,'=futurelattice'

from __future__ import division
from scipy.ndimage import measurements as sci
import numpy as np
import functions
import functionsdaan
#import ani
import matplotlib.pyplot as plt


size =20 #Width of the lattice
init_M = 0.99 #initial magnetization
lattice = functions.initialize(size,init_M) #initialize the girl, M sets magnetization
while np.average(lattice) != init_M:
    lattice =functions.initialize(size,init_M)
T = 0.4  #temperature, entropy vs potential energy
invT = 1./T
beta = invT
runtime = 5000 #Total time of the simulation
B = 0.1 #The favorite orientation
Trange = np.linspace(0.001,1,5)

#Initializing arrays for clusterflippin'
blacklist = np.zeros((len(lattice),len(lattice[0])))
futurelattice =np.zeros((len(lattice),len(lattice[0])))
xval = 5 #np.random.randint(0,len(lattice)+1)
yval = 5 #np.random.randint(0,len(lattice)+1)
blacklist[yval][xval] = 1
lattice[yval][xval] = -1*lattice[yval][xval] #the module about accepting or rejecting the first randomly picked spin must be built still
futurelattice[yval][xval] = lattice[yval][xval]
comp_xval = xval + 1
comp_yval = yval
blacklist[comp_yval][comp_xval] = 1
cnt = 0
futurelattice = functionsdaan.wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)

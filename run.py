# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:39:10 2015

@author: Willem
"""


from __future__ import division
from scipy.ndimage import measurements as sci
import numpy as np
import functions
#import ani
import matplotlib.pyplot as plt

#Investigate the behaviour of an Ising lattice using a Markov chain;
#flipping random spins and adopting the new situation with a boltzman probality.

#Set the temperature, size of the lattice, the runtime of the simulation,
#and the initial magnetization, and out comes
#the magnetization after the runtime!
plt.close("all")

size =20 #Width of the lattice
init_M = 0.99 #initial magnetization
init_l = functions.initialize(size,init_M) #initialize the girl, M sets magnetization
T = 0.5  #temperature, entropy vs potential energy
invT = 1/T
runtime = 5000 #Total time of the simulation
B = 0.1 #The favorite orientation
Trange = np.linspace(0.001,1,5)
critval = 0.15
flipnum = 2
test=-1*np.ones((9,9))
test[::2,0]= 1
test[::2,-1]= 1 #BOTLEFT
test[0,::2]=1
test[-1,::2]=1

switch = 3
if switch == 1:
    lattice = init_l
elif switch == 2:
    lattice = functions.simulate(init_l,invT,B,50000)
    print "ooo"
elif switch == 3:
  lattice = init_l
  while np.average(lattice) != init_M:
    lattice =functions.initialize(size,init_M)
  Marray = functions.crittemp(lattice,Trange,critval,flipnum,runtime)
else:
    lattice = test
    
print Marray
    
plt.plot(Marray,range(len(Marray)))
plt.show()
    
"""  
plt.figure(1)
plt.subplot(121)
plt.imshow(lattice,cmap="Greys",interpolation="nearest")
plt.colorbar()
plt.subplot(122)
plt.imshow(functions.clusters(lattice)[0],cmap="Set2",interpolation="nearest")
plt.colorbar()
plt.show()




#run the simulation for runtime steps, at temp T and B-field B
print "B: ",B
print "T: ",T
lattice = functions.simulate(init_Lattice,invT,B,runtime)
M = np.average(lattice) #calc et print magnetization
print M
plt.imshow(lattice, cmap = "Greys",interpolation = "nearest")
"""


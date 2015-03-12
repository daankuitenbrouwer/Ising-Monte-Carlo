# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:39:10 2015

@author: Willem
"""


from __future__ import division
from scipy.ndimage import measurements as sci
import numpy as np
import functions
import matplotlib.pyplot as plt

#Investigate the behaviour of an Ising lattice using a Markov chain;
#flipping random spins and adopting the new situation with a boltzman probality.

#Set the temperature, size of the lattice, the runtime of the simulation,
#and the initial magnetization, and out comes
#the magnetization after the runtime!
plt.close("all")

size = 100 #Width of the lattice
init_M = 0 #initial magnetization
lattice = functions.initialize(size,init_M) #initialize the girl, M sets magnetization
T = 100  #temperature, entropy vs potential energy
invT = 1/T
runtime = 10000 #Total time of the simulation
B = 0 #The favorite orientation
runs = 10


for i in range(runs):
    init_M = 2*i/runs - 1
    lattice = functions.initialize(size,init_M)
    
    print i
    M = functions.loop(lattice,invT,B,runtime)
    plt.plot(M)
    


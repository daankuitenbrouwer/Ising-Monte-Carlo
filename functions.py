# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:38:41 2015

@author: Willem
"""

from __future__ import division
import numpy as np

#initializes a random lattice.
def initialize(size,init_M):
    bin_lattice = np.random.random((size,size))+0.5*init_M
    lattice = 2*np.round(bin_lattice)-1
    return lattice

#Takes a lattice and return an array containg the energies of all the spins.
def global_energy(lattice,B):
    l_shifted = np.roll(lattice,1,axis=1)
    r_shifted = np.roll(lattice,-1,axis=1)
    d_shifted = np.roll(lattice,1,axis=0)
    u_shifted = np.roll(lattice,-1,axis=0)    
    
    #minus sign to make it antiferromagnetic 
    energy_array = -1*(lattice*l_shifted + lattice*r_shifted + \
    lattice*u_shifted + lattice*d_shifted) + B*lattice
    
    return energy_array

#Flips a random spin at coord. x,y, calculates the caused difference in energy.  
def flip_and_dE(lattice,B):
    width = len(lattice)
    x = np.random.randint(0,width)
    y = np.random.randint(0,width)
    old_Energy = calc_energy_local(lattice,B,x,y)
    temp_lattice = np.empty_like (lattice) #deep copy the old lattice
    temp_lattice[:] = lattice    
    
    temp_lattice[x,y] = -1*temp_lattice[x,y] #flip the spin
    
    new_Energy = calc_energy_local(temp_lattice,B,x,y)
    dE = new_Energy - old_Energy
    return temp_lattice,dE,x,y

#Calculates the energy of a central spin and its 4 adjacent neighbours. 
def calc_energy_local(lattice,B,x,y):
    
    if (x == len(lattice)-1):
        d = lattice[0,y]
    else:
        d = lattice[x+1,y]
        
    if (y == len(lattice)-1):
        r = lattice[x,0]    
    else:
        r = lattice[x,y+1]
    
    u = lattice[x-1,y]
    l = lattice[x,y-1]
    
    s = lattice[x,y]
    
    #the energy of all pairs + the energy of the single spin in the B-field
    loc_energy = -2*(s * (l + r + u + d)) - s*B
    
    #debugging
    #loc_energy = -s*B
    return loc_energy

#This script takes two lattices, the temperature and the difference in Energy
#Between the lattices. A chance is calculated (exp(-dE/T)) with which 
#the new lattice is adopted.
def take_or_refuse(invT,old_lattice,new_lattice,dE):
    P = np.exp(-dE*invT)
    r_num = np.random.random(1)
    #print "r_num: ",r_num
    is_chosen = 0
    if P>= 1:
        choice = new_lattice
        is_chosen = 1
    elif r_num< P:
        choice = new_lattice
        is_chosen = 1
    else:
        choice = old_lattice
        is_chosen = 0
    #print "Is chosen:", is_chosen, "with P: ",P
    return choice,P,is_chosen
 
#Lets the simulation run, for a duration of runtime. Applies a B field
#during the startup.   
def loop(lattice,invT,B,runtime):
    M = np.zeros(runtime) #magnetization
    for i in range(runtime):
        
        M[i] = np.average(lattice)
        new_lattice,dE,x,y = flip_and_dE(lattice,B)    
        
        lattice,P,is_chosen = take_or_refuse(invT,lattice,new_lattice,dE)        
         
    return M
        
    

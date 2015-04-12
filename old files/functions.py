# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 15:38:41 2015

@author: Willem
"""

from __future__ import division
import numpy as np
from pylab import *
from scipy.ndimage import measurements
import matplotlib.pyplot as plt
import random

#initializes a random lattice.
def initialize(size,init_M):
    pre_lattice = np.random.random((size,size))
    magprelattice = np.average(pre_lattice)
    print magprelattice, '=magprelattice'
    bin_lattice = pre_lattice - (magprelattice - 0.5) + 0.5*init_M
    lattice = 2*np.round(bin_lattice)-1
    return lattice.astype(int)

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
def simulate(lattice,invT,B,runtime):
    M = np.zeros(runtime) #magnetization
    for i in range(runtime):
        
        M[i] = np.average(lattice)
        new_lattice,dE,x,y = flip_and_dE(lattice,B)    
        
        lattice,P,is_chosen = take_or_refuse(invT,lattice,new_lattice,dE)        
         
    return lattice

def step(lattice,invT,B):
    new_lattice,dE,x,y = flip_and_dE(lattice,B)    
        
    lattice,P,is_chosen = take_or_refuse(invT,lattice,new_lattice,dE)        
         
    return lattice


#return naive clusters areas
def cluster_calc(lattice):
     lw, num = measurements.label(0.5*(lattice+1))
     area = measurements.sum(lattice, lw, index=arange(lw.max() + 1))
     areaImg = area[lw]
     return area,areaImg,lw

#Takes labels array and returns same structure with lowest ints instead.
def improve_label(labels):
    mapper = np.arange(0,np.max(labels)+1,1)
    key = np.unique(labels)
    for i in range(len(key)):
        mapper[key[i]]=i
    new_label = mapper[labels]
    return new_label
    
#The scipy package returns a an array with labeled clusters, and the labels,
#but doesnt take into account the periodic boundary conditions.
#This script takes into account the periodic bcs, and looks only at up-spin!
#It adds all labels that are in fact on the same cluster and sets them equal
#Then uses these labels to set the new clusters array
def sci_pbc(lattice):
    c_lattice = measurements.label(0.5*(lattice+1))[0]
    labels = np.arange(0,measurements.label(0.5*(lattice+1))[1]+1,1)
    #First check adjacent corners
    labels = test_update(lattice, c_lattice,labels,0,0,-1,0)
    labels = test_update(lattice, c_lattice,labels,-1,0,-1,-1)
    labels = test_update(lattice, c_lattice,labels,-1,-1,0,0-1)
    labels = test_update(lattice, c_lattice,labels,0,-1,0,0)
    #Check the vertical sides
    for i in range(len(lattice)-2):
        labels = test_update(lattice, c_lattice,labels,0,i+1,-1,i+1) 
    #Check the horizontal sides
    for i in range(len(lattice)-2):
        labels = test_update(lattice, c_lattice,labels,i+1,0,i+1,-1)
    #update labels to take into account periodic bcs
    labels = improve_label(labels)
    new_clusters = labels[c_lattice]
    return new_clusters,labels.max()

#Returns array containing labeled clusters and highest label value.
def cluster_array(lattice):
    #mask1 = (0.5*(lattice+1)).astype(int)
    cluster1,num1 = sci_pbc(lattice)    
    #mask2 = (-0.5*(lattice-1)).astype(int)
    cluster2,num2= sci_pbc(-1*lattice)
    new_labels = np.arange(num1-1,num1+num2,1)
    new_labels[0] = 0
    return cluster1 + new_labels[cluster2] -1,num1+num2

#Labels label clusters. Cluster 1 has label 1 , cluster 2 has label2 etc.
#Takes a labels array and two labels and sets the labels equal to each other.
def label_update2(labels,label1,label2):
    mapper = np.arange(0,len(labels),1)
    #set values in map equal to lowest of label1,label2
    mapper[label1]=np.amin(np.array(labels[label1],labels[label2]))
    mapper[label2]=np.amin(np.array(labels[label1],labels[label2]))
    new_label = mapper[labels]
    return new_label    
    
#Takes two coordinates,a labels array and a 01 spin lattice 
#checks if the lattice has same spin at those coords,
#if so it sets the labels of the clusters at those coords
#equal.
def test_update(lattice,c_lattice,labels,x1,y1,x2,y2):
    if(lattice[x1,y1]==lattice[x2,y2]):
        return label_update2(labels,c_lattice[x1,y1],c_lattice[x2,y2])    
    else: 
        return labels
    
def colorful_cluster_visuals(lattice):
    lw, num = measurements.label(0.5*(lattice+1))
    print lw
    print num
    b = arange(lw.max() + 1) # create an array of values from 0 to lw.max() + 1
    print b  
    shuffle(b) # shuffle this array
    print b
    shuffledLw = b[lw] # replace all values with values from b
    print shuffledLw
    #imshow(shuffledLw, origin='lower', interpolation='nearest')
     
def cluster_hist(area,bins):
    cluster_hist = np.histogram(area,bins)    
    plt.plot(np.linspace(np.min(area),np.max(area),bins),cluster_hist[0])
    return cluster_hist[0],np.min(area),np.max(area)
    
   
def clusterflipper(lattice):  
  clusterarray,highlabel = cluster_array(lattice)
  labelnum = np.random.randint(0,highlabel+1)
  print clusterarray, '=clusterarray',labelnum,'=labelnum'
  fliploc = np.where(clusterarray == labelnum)
  print fliploc, '=fliploc'
  for coord in fliploc:
    lattice[coord] = -1*lattice[coord]
  return lattice
    
def crittemp(lattice,Trange,critval,flipnum,runtime):
  Marray = np.zeros((len(Trange),flipnum))
  
  for T in Trange:
    print np.average(lattice), 'initialm'
    invT = (1./T)
    print T,'=T'
    for flp in range(flipnum):
      print np.average(lattice), '=m'
      formerlattice = np.empty_like(lattice)
      formerlattice[:] = lattice
      lattice = simulate(lattice,invT,0,runtime)
      lattice = clusterflipper(lattice)
      plt.figure(1)
      plt.subplot(121)
      plt.imshow(lattice,cmap="Greys",interpolation="nearest")
      plt.colorbar()
      plt.subplot(122)
      plt.imshow(formerlattice,cmap="Greys",interpolation="nearest")
      plt.colorbar()
      plt.show()
      print formerlattice, 'formerlattice',lattice,'lattice'
      dE = np.sum(global_energy(formerlattice,0)) - np.sum(global_energy(lattice,0))
      lattice,P,is_chosen = take_or_refuse(invT,formerlattice,lattice,dE)
      Marray[flp] = np.average(lattice)
      print Marray[T][flp],'=M',T,'=T'
  return Marray
  

  
  




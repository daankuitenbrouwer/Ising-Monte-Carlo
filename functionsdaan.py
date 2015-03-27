from __future__ import division
import numpy as np
from pylab import *
from scipy.ndimage import measurements
import matplotlib.pyplot as plt
import random
from sys import exit

def initializelat(size,init_M):
  pre_lattice = np.random.random((size,size))
  magprelattice = np.average(pre_lattice)
  #print magprelattice, '=magprelattice'
  bin_lattice = pre_lattice - (magprelattice - 0.5) + 0.5*init_M
  lattice = 2*np.round(bin_lattice)-1
  return lattice.astype(int)
  
def initializeelse(size):
  blacklist = np.zeros((size,size))
  futurelattice = np.zeros((size,size))
  return blacklist,futurelattice

def rand_init_pos(blacklist,lattice):
  #print blacklist
  available_lst = np.where(blacklist == 0)
  #print available_lst,'=avabl',len(available_lst),len(available_lst[0])
  rnl = np.linspace(0,len(available_lst[0])-1,len(available_lst[0]))
  picked_val = random.choice(rnl)
  #print picked_val,'=pval',np.max(rnl)
  xval = available_lst[1][picked_val]
  yval = available_lst[0][picked_val]
  #print xval,yval,'=x,y from rand_init_pos'
  return xval,yval

def take_reject_neighbour(lattice,xval,yval,comp_xval,comp_yval,blacklist,beta):
#  if lattice[yval][xval] == lattice[comp_yval][comp_xval]:
#   acceptflip = 1 # This means that the flip is not accepted
  #R = random.uniform(-1,1)
  P = 1-np.exp(np.min([0,2*beta*lattice[yval][xval]*lattice[comp_yval][comp_xval]]))
  compval = random.uniform(0,1)
  #print P,compval,beta,'=P,compval,beta,',lattice[yval][xval]*lattice[comp_yval][comp_xval],'sisj'
  if P > compval:
    acceptflip = -1
  else:
    acceptflip = 1
  if blacklist[comp_yval][comp_xval] != 0:
    #print 'here',blacklist[comp_yval][comp_xval],'=bl'
    acceptflip = 1
  return acceptflip
  
def periodic_boundary(lattice,xdir,ydir,xval,yval,comp_xval,comp_yval):
  dir_mat = np.array([xdir,ydir])
  coord_matrix = np.array([xval,yval,comp_xval,comp_yval])
  SW_outside = np.where(coord_matrix == len(lattice))
  NE_outside = np.where(coord_matrix == -1)
  coord_matrix[SW_outside] = 0
  coord_matrix[NE_outside] = len(lattice) - 1
  positive_movement = np.where(dir_mat == -(len(lattice) - 1))
  negative_movement = np.where(dir_mat == len(lattice) - 1)
  dir_mat[positive_movement] = 1
  dir_mat[negative_movement] = -1
  
  #print dir_mat,'=dir_mat', coord_matrix,'=coord_matrix'
  xdir = dir_mat[0]
  ydir = dir_mat[1]
  xval = coord_matrix[0]
  yval = coord_matrix[1]
  comp_xval = coord_matrix[2]
  comp_yval = coord_matrix[3]
  #print xval,yval,comp_xval,comp_yval,'=xval,yval,comp_xval,comp_yval should be coord_matrix'
  return xdir,ydir,xval,yval,comp_xval,comp_yval
  
def blackboxchecker(xdir,ydir,xval,yval,comp_xval,comp_yval,blacklist):
  #print 'BBC'
  newseed = 0
  if blacklist[comp_yval][comp_xval] != 0:
    comp_xval = xval + xdir
    comp_yval = yval + ydir
    xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(blacklist,xdir,ydir,xval,yval,comp_xval,comp_yval)
    if blacklist[comp_yval][comp_xval] != 0:
      comp_xval = xval + ydir
      comp_yval = yval - xdir
      xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(blacklist,xdir,ydir,xval,yval,comp_xval,comp_yval)
      if blacklist[comp_yval][comp_xval] != 0:
        newseed = 1
  return xval,yval,comp_xval,comp_yval,ydir,xdir,newseed    
  
def position(acceptflip,blacklist,xval,yval,comp_xval,comp_yval):
  xdir = comp_xval - xval
  ydir = comp_yval - yval
  xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(blacklist,xdir,ydir,xval,yval,comp_xval,comp_yval)
  if acceptflip == -1:
    xval = comp_xval
    yval = comp_yval
    xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(blacklist,xdir,ydir,xval,yval,comp_xval,comp_yval)
    comp_xval = xval - ydir
    comp_yval = yval + xdir
    xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(blacklist,xdir,ydir,xval,yval,comp_xval,comp_yval)
  else:
    comp_xval = xval + ydir
    comp_yval = yval - xdir
    xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(blacklist,xdir,ydir,xval,yval,comp_xval,comp_yval)
    ydir = -ydir
    xdir = -xdir
    #xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(blacklist,xdir,ydir,xval,yval,comp_xval,comp_yval)
  #print xval,comp_xval,yval,comp_yval,'pos from pbc'
  #print xdir,ydir,'=dirs'
  xval,yval,comp_xval,comp_yval,ydir,xdir,newseed = blackboxchecker(xdir,ydir,xval,yval,comp_xval,comp_yval,blacklist)
  #print xdir,ydir,'=dirs'
  return xval,yval,comp_xval,comp_yval,ydir,xdir,newseed
    
def wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta):
  cnt = cnt +1
  if cnt > 994:
    print cnt, 'for cmp check'
    exit(0)
  acceptflip = take_reject_neighbour(lattice,xval,yval,comp_xval,comp_yval,blacklist,beta)
  #print acceptflip,'=acceptflip',lattice[comp_yval][comp_xval],lattice[yval][xval],'=lattice[comp_xval],lattice[xval]',comp_xval,'comp_xval'
  blacklist[comp_yval][comp_xval] = 1
  blacklist[yval][xval] = blacklist[yval][xval] + 1
  lattice[comp_yval][comp_xval] = acceptflip*lattice[comp_yval][comp_xval]
  futurelattice[comp_yval][comp_xval] = lattice[comp_yval][comp_xval]
  xval,yval,comp_xval,comp_yval,ydir,xdir,newseed = position(acceptflip,blacklist,xval,yval,comp_xval,comp_yval)
  #print len(np.where(blacklist == 0)[0]), '#not attended sites'
  if np.min(blacklist) == 0:
    '''
    #print xval,'=xval',comp_xval,'=comp_xval',yval,'=yval',comp_yval,'=comp_yval'
    cnt = cnt + 1
    plt.figure(1)
    plt.subplot(131)
    plt.imshow(lattice,cmap="Greys",interpolation="nearest")
    #plt.colorbar()
    plt.subplot(132)
    plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
    #plt.colorbar()
    plt.subplot(133)
    plt.imshow(blacklist,cmap="Set2",interpolation="nearest")
    #plt.colorbar()
    if cnt%3 == 0:
      plt.show()
    #newseed = 1
     '''  
    if newseed == 1:
      #print np.min(blacklist),'=min blacklist'
      while newseed == 1 and blacklist[yval][xval] != 0:
        xval,yval = rand_init_pos(blacklist,lattice) # New chosen flip is unconditionally accepted, need single flip method to decide upon accepting or rejecting?
      
      blacklist[yval][xval] = 1
      
      
      lattice[yval][xval] = -1*lattice[yval][xval] #the module about accepting or rejecting the first randomly picked spin must be built still
      futurelattice[yval][xval] = lattice[yval][xval]
      comp_xval = xval + 1
      comp_yval = yval
      xdir = comp_xval - xval
      ydir = comp_yval - yval
      xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(lattice,xdir,ydir,xval,yval,comp_xval,comp_yval)
      xval,yval,comp_xval,comp_yval,ydir,xdir,newseed = blackboxchecker(xdir,ydir,xval,yval,comp_xval,comp_yval,blacklist)
      
    wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)
  return futurelattice
  
def crittemp(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,Trange,size):
  beta = 1./Trange[0]
  cnt = 0
  futurelattice = wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)
  magnetization = np.zeros(len(Trange))
  
  for ind,T in enumerate(Trange):
    beta = 1./T
  #plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
  #print 'fut'
  #plt.show()
    lattice = np.empty_like(futurelattice)
    lattice[:] = futurelattice
  #plt.imshow(lattice,cmap="Greys",interpolation="nearest")
  #print 'lat'
  #plt.show()
    blacklist,futurelattice = initializeelse(size)
    xval,yval = rand_init_pos(blacklist,lattice) 
    blacklist[yval][xval] = 1 # I flip the first spin unconditionally
    lattice[yval][xval] = -1*lattice[yval][xval] #the module about accepting or rejecting the first randomly picked spin must be built still
    futurelattice[yval][xval] = lattice[yval][xval]
    comp_xval = xval + 1
    comp_yval = yval
    xdir = comp_xval - xval
    ydir = comp_yval - yval
    xdir,ydir,xval,yval,comp_xval,comp_yval = periodic_boundary(lattice,xdir,ydir,xval,yval,comp_xval,comp_yval)
    cnt = 0
  #print lattice[yval][xval],lattice[comp_xval][comp_yval],'latxy,latcompxy'
    futurelattice = wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)
    magnetization[ind] = np.average(futurelattice)
    futurelattice[:] = -1*futurelattice[:]
  #print magnetization[ind],'magnetization[ind]'
  #plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
  #plt.show()
#  print magnetization
  
  return magnetization

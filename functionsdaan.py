from __future__ import division
import numpy as np
from pylab import *
from scipy.ndimage import measurements
import matplotlib.pyplot as plt
import random


def take_reject_neighbour(lattice,xval,yval,comp_xval,comp_yval,blacklist,beta):
  R = random.uniform(-1,1)
  P = 1-np.exp(np.min([0,2*beta*R*R*lattice[yval][xval]*lattice[comp_yval][comp_xval]]))
  compval = random.uniform(0,1)
  if P > compval:
    acceptflip = -1
  else:
    acceptflip = 1
  if blacklist[comp_yval][comp_xval] != 0:
    acceptflip = 1
  return acceptflip
  
def blackboxchecker(xdir,ydir,xval,yval,comp_xval,comp_yval,blacklist):
  if blacklist[comp_yval][comp_xval] == 1:
    comp_xval = xval + xdir
    comp_yval = yval + ydir
    if blacklist[comp_yval][comp_xval] == 1:
      comp_xval = xval + ydir
      comp_yval = yval - xdir
      if blacklist[comp_yval][comp_xval] == 1:
        xval = xval - xdir
        yval = yval - ydir
        comp_xval = xval + ydir
        comp_yval = yval - xdir
        ydir = -ydir
        xdir = -xdir
  return xval,yval,comp_xval,comp_yval,ydir,xdir    
  
def position(acceptflip,blacklist,xval,yval,comp_xval,comp_yval):
  xdir = comp_xval - xval
  ydir = comp_yval - yval
  if acceptflip == -1:
    xval = comp_xval
    yval = comp_yval
    comp_xval = xval - ydir
    comp_yval = yval + xdir
  else:
    comp_xval = xval + ydir
    comp_yval = yval - xdir
    ydir = -ydir
    xdir = -xdir
  xval,yval,comp_xval,comp_yval,ydir,xdir = blackboxchecker(xdir,ydir,xval,yval,comp_xval,comp_yval,blacklist)
  acceptflip
  return xval,yval,comp_xval,comp_yval,ydir,xdir
    
def wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta):
  acceptflip = take_reject_neighbour(lattice,xval,yval,comp_xval,comp_yval,blacklist,beta)
  blacklist[comp_yval][comp_xval] = 1
  blacklist[yval][xval] = blacklist[yval][xval] + 1
  lattice[comp_yval][comp_xval] = acceptflip*lattice[comp_yval][comp_xval]
  futurelattice[comp_yval][comp_xval] = lattice[comp_yval][comp_xval]
  xval,yval,comp_xval,comp_yval,ydir,xdir = position(acceptflip,blacklist,xval,yval,comp_xval,comp_yval)
  if np.min(blacklist) == 0:
    plt.figure(1)
    plt.subplot(131)
    plt.imshow(lattice,cmap="Greys",interpolation="nearest")
    plt.colorbar()
    plt.subplot(132)
    plt.imshow(futurelattice,cmap="Greys",interpolation="nearest")
    plt.colorbar()
    plt.subplot(133)
    plt.imshow(blacklist,cmap="Set2",interpolation="nearest")
    plt.colorbar()
    plt.show()
    wolff(lattice,blacklist,futurelattice,xval,yval,comp_xval,comp_yval,cnt,beta)
  return futurelattice

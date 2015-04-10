import numpy as np
import matplotlib.pyplot as plt
import random 

def InitializeLattice(size,init_M):
  pre_lattice = np.random.random((size,size))
  magprelattice = np.average(pre_lattice)
  bin_lattice = pre_lattice - (magprelattice - 0.5) + 0.5*init_M
  lattice = 2*np.round(bin_lattice)-1
  while np.average(lattice.astype(int)) != init_M:
    lattice,BlackList = InitializeLattice(size,init_M)
  BlackList = np.ones((size,size),dtype=bool)
  return lattice.astype(int),BlackList
  
def Positions_Neighbours(y,x,lattice):
  up = (y - 1)%(len(lattice))
  down = (y + 1)%(len(lattice))
  left = (x - 1)%(len(lattice))
  right = (x + 1)%(len(lattice))
  return up,down,left,right

def GrowCluster(y,x,lattice,BlackList,J,cnt):
  up,down,left,right = Positions_Neighbours(y,x,lattice)
  BlackList[y][x] = False #In the beginning everything should be 'True'
  cnt = cnt + 1 # counter can be used for visualisation, see plotting scheme below
  if BlackList[up][x]:
    lattice, BlackList,cnt = TryAdd(up,x,y,x,lattice,BlackList,J,cnt)
  if BlackList[y][right]:
    lattice, BlackList,cnt = TryAdd(y,right,y,x,lattice,BlackList,J,cnt)
  if BlackList[down][x]:
    lattice, BlackList,cnt = TryAdd(down,x,y,x,lattice,BlackList,J,cnt)
  if BlackList[y][left]:
    lattice, BlackList,cnt = TryAdd(y,left,y,x,lattice,BlackList,J,cnt)
  else:
    pass
  return lattice,BlackList,cnt
  

def TryAdd(TryAdd_y,TryAdd_x,y,x,lattice,BlackList,J,cnt):
  spin_val_difference = lattice[y][x]*lattice[TryAdd_y][TryAdd_x]
  if spin_val_difference == 1:
    print spin_val_difference,'svd'
  P = 1. - np.exp(np.min([0.,2.*J*spin_val_difference]))
  compval = random.uniform(0,1)
  #print P, compval,(1./J),'P,compval,T'
  if P > compval:
    lattice[TryAdd_y][TryAdd_x] = -1*lattice[TryAdd_y][TryAdd_x]
    GrowCluster(TryAdd_y,TryAdd_x,lattice,BlackList,J,cnt)
  return lattice, BlackList,cnt
  
def Cluster_Creator(lattice,J,Niter):
  BlackList = np.ones((len(lattice),len(lattice)),dtype=bool)
  cnt = 0
  for i in range(Niter):
    y = random.choice(BlackList[0])
    x = random.choice(BlackList[0])
    lattice[y][x] = -1*lattice[y][x]
    lattice,BlackList,cnt = GrowCluster(y,x,lattice,BlackList,J,cnt)
  return lattice
  
def CritTemp(lattice,BlackList,Trange,Niter):
  Magnetization_func_of_T = np.zeros(len(Trange))
  for ind,T in enumerate(Trange):
    J = 1./T
    lattice = Cluster_Creator(lattice,J,Niter)
    Magnetization_func_of_T[ind] = np.average(lattice)
  return Magnetization_func_of_T
  
'''
    plt.figure(1)
    plt.subplot(121)
    plt.imshow(lattice,cmap="Greys",interpolation="nearest")
    plt.subplot(122)
    plt.imshow(BlackList.astype(int),cmap="Greys",interpolation="nearest")
    plt.show()
     '''

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
  
def rand_init_pos(BlackList):
  available_lst = np.where(BlackList)
  available_indices = np.linspace(0,len(available_lst[0])-1,len(available_lst[0]))
  picked_index = random.choice(available_indices)
  x = available_lst[1][picked_index]
  y = available_lst[0][picked_index]
  return y,x
  
def Positions_Neighbours(y,x,lattice):
  up = (y - 1)%(len(lattice))
  down = (y + 1)%(len(lattice))
  left = (x - 1)%(len(lattice))
  right = (x + 1)%(len(lattice))
  return up,down,left,right

def GrowCluster(y,x,lattice,BlackList,J,cnt):
  up,down,left,right = Positions_Neighbours(y,x,lattice)
  BlackList[y][x] = False #In the beginning everything should be 'True'
  cnt = cnt + 1
  '''
  if cnt%6 == 0:
    plt.figure(1)
    plt.subplot(121)
    plt.imshow(lattice,cmap="Greys",interpolation="nearest")
    plt.subplot(122)
    plt.imshow(BlackList.astype(int),cmap="Greys",interpolation="nearest")
    plt.show()
   '''
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
  BlackList[TryAdd_y][TryAdd_x] = False
  spin_val_difference = lattice[y][x]*lattice[TryAdd_y][TryAdd_x]
  P = 1. - np.exp(np.min([0.,2.*J*spin_val_difference]))
  compval = random.uniform(0,1)
  #print P, compval,(1./J),'P,compval,T'
  if P > compval:
    lattice[TryAdd_y][TryAdd_x] = -1*lattice[TryAdd_y][TryAdd_x]
    GrowCluster(TryAdd_y,TryAdd_x,lattice,BlackList,J,cnt)
  return lattice, BlackList,cnt
  
def accept_refuse_new_flip(y,x,lattice,J):
  #up,down,left,right = Positions_Neighbours(y,x,lattice)
  #dE = 2*-lattice[y][x]*(lattice[up][x] + lattice[down][x] + lattice[y][left] + lattice[y][right])
  P = 1#np.exp(2*J*dE)
  compval = random.uniform(0,1)
  #print P, compval 
  acceptflip = 0
  if P > compval:
    lattice[y][x] = -1*lattice[y][x]
    acceptflip = 1
  return lattice,acceptflip
  
def Cluster_Creator(lattice,J):
  BlackList = np.ones((len(lattice),len(lattice)),dtype=bool)
  number_non_visited = len(np.where(BlackList)[0])
  cnt = 0
  while number_non_visited != 0:
    y,x = rand_init_pos(BlackList)
    
    BlackList[y][x] = False
    lattice, acceptflip =accept_refuse_new_flip(y,x,lattice,J)
    '''
    plt.figure(1)
    plt.subplot(121)
    plt.imshow(lattice,cmap="Greys",interpolation="nearest")
    plt.subplot(122)
    plt.imshow(BlackList.astype(int),cmap="Greys",interpolation="nearest")
    plt.show()
     '''
    if acceptflip == 0:
      number_non_visited = len(np.where(BlackList)[0])
      continue
    else:
      lattice, BlackList,cnt = GrowCluster(y,x,lattice,BlackList,J,cnt)
      number_non_visited = len(np.where(BlackList)[0])
  return lattice
  
def CritTemp(lattice,BlackList,Trange):
  Magnetization_func_of_T = np.zeros(len(Trange))
  for ind,T in enumerate(Trange):
    BlackList = np.ones((len(BlackList),len(BlackList)),dtype=bool)
    J = 1./T
    lattice = Cluster_Creator(lattice,J)
    Magnetization_func_of_T[ind] = np.average(lattice)
    lattice = -1*lattice
  return Magnetization_func_of_T
  


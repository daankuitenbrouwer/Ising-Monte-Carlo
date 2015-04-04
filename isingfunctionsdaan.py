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

def GrowCluster(y,x,lattice,BlackList,J,Rec_cnt,maxrecursion,cnt):
  up,down,left,right = Positions_Neighbours(y,x,lattice)
  BlackList[y][x] = False #In the beginning everything should be 'True'
  Rec_cnt = Rec_cnt + 1 #Recursion counter
  #print Rec_cnt,maxrecursion,'Rec_cnt,maxrecursion',cnt
  cnt = cnt + 1
  '''
  if cnt%(10) == 0:
    plt.figure(1)
    plt.subplot(121)
    plt.imshow(lattice,cmap="Greys",interpolation="nearest")
    plt.subplot(122)
    plt.imshow(BlackList.astype(int),cmap="Greys",interpolation="nearest")
    plt.show()
   '''
  if Rec_cnt < maxrecursion:
    if BlackList[up][x]:
      lattice, BlackList = TryAdd(up,x,y,x,lattice,BlackList,J,Rec_cnt,maxrecursion,cnt)
    if BlackList[y][right]:
      lattice, BlackList = TryAdd(y,right,y,x,lattice,BlackList,J,Rec_cnt,maxrecursion,cnt)
    if BlackList[down][x]:
      lattice, BlackList = TryAdd(down,x,y,x,lattice,BlackList,J,Rec_cnt,maxrecursion,cnt)
    if BlackList[y][left]:
      lattice, BlackList = TryAdd(y,left,y,x,lattice,BlackList,J,Rec_cnt,maxrecursion,cnt)
  else:
    pass
  return lattice,BlackList, Rec_cnt
  

def TryAdd(TryAdd_y,TryAdd_x,y,x,lattice,BlackList,J,Rec_cnt,maxrecursion,cnt):
  BlackList[TryAdd_y][TryAdd_x] = False
  spin_val_difference = lattice[y][x]*lattice[TryAdd_y][TryAdd_x]
  P = 1-np.exp(np.min([0,2*J*spin_val_difference]))
  compval = random.uniform(0,1)
  if P > compval:
    lattice[TryAdd_y][TryAdd_x] = -1*lattice[TryAdd_y][TryAdd_x]
    GrowCluster(TryAdd_y,TryAdd_x,lattice,BlackList,J,Rec_cnt,maxrecursion,cnt)
  return lattice, BlackList
  
def accept_refuse_new_flip(y,x,lattice,J):
  up,down,left,right = Positions_Neighbours(y,x,lattice)
  dE = -2*-lattice[y][x]*(lattice[up][x] + lattice[down][x] + lattice[y][left] + lattice[y][right])
  P = np.exp(-2*dE*J)
  compval = random.uniform(0,1) 
  acceptflip = 0
  if P > compval:
    lattice[y][x] = -1*lattice[y][x]
    acceptflip = 1
  return lattice,acceptflip
  
def Cluster_Creator(lattice,J,maxrecursion):
  BlackList = np.ones((len(lattice),len(lattice)),dtype=bool)
  number_non_visited = len(np.where(BlackList)[0])
  while number_non_visited != 0:
    y,x = rand_init_pos(BlackList)
    BlackList[y][x] = False
    lattice, acceptflip =accept_refuse_new_flip(y,x,lattice,J)
    if acceptflip == 0:
      number_non_visited = len(np.where(BlackList)[0])
      continue
    else:
      Rec_cnt = 0
      cnt = 0
      lattice, BlackList, Rec_cnt = GrowCluster(y,x,lattice,BlackList,J,Rec_cnt,maxrecursion,cnt)
      number_non_visited = len(np.where(BlackList)[0])
       
  return lattice
  
def CritTemp(lattice,BlackList,Trange,maxrecursion):
  Magnetization_func_of_T = np.zeros(len(Trange))
  for ind,T in enumerate(Trange):
    BlackList = np.ones((len(BlackList),len(BlackList)),dtype=bool)
    J = 1./T
    lattice = Cluster_Creator(lattice,J,maxrecursion)
    Magnetization_func_of_T[ind] = np.average(lattice)
    lattice = -1*lattice
  return Magnetization_func_of_T
  


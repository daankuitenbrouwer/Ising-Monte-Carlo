import numpy as np
import random as rn
import scipy


def init_Lattice(N,Nstep):
  Lattice = 2*np.random.randint(2,size=(N,N))-1 #how to exclude the zeros here?
  #Lattice = np.zeros((N,N)) + 1
  marray = np.zeros(Nstep, dtype=float)
  marray[0] = np.sum(Lattice)
  parray = np.zeros(Nstep, dtype=float)
  dearray = np.zeros(Nstep, dtype=float)
  return Lattice, marray, parray, dearray
  
def initial_energy(Lattice,N):
#replacing the matrix vertically
  upside = Lattice[0,:]
  downside = Lattice[-1,:]
  intermedUP = scipy.delete(Lattice,0,0) 
  LatticeU = np.concatenate((intermedUP,[upside]),axis=0)
  intermedDO = scipy.delete(Lattice,-1,0) 
  LatticeDO = np.concatenate(([downside],intermedDO),axis=0)

#replacing the matrix horizontally
  rotatedlat = np.array(zip(*Lattice[::-1]))
  leftside = rotatedlat[0,:]
  rightside = rotatedlat[-1,:]
  intermedLE = scipy.delete(rotatedlat,0,0) 
  LatticeLE = np.concatenate((intermedLE,[leftside]),axis=0)
  intermedRI = scipy.delete(rotatedlat,-1,0) 
  LatticeRI = np.concatenate(([rightside],intermedRI),axis=0)

#subtracting and working out: SUM (Si*Sj)
  upsubtr = Lattice - LatticeU
  dosubtr = Lattice - LatticeDO
  lesubtr = rotatedlat - LatticeLE
  risubtr = rotatedlat - LatticeRI


  upcontribution = -abs(upsubtr)+1
  docontribution = -abs(dosubtr)+1
  lecontribution = -abs(lesubtr)+1
  ricontribution = -abs(risubtr)+1
  totcontribution = np.sum(upcontribution + docontribution + lecontribution + ricontribution)
  print totcontribution/(N**2.), '=initialmagnetization'
  return totcontribution

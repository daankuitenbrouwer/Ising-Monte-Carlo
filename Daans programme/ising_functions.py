import numpy as np

#accepting or declining the new lattice
def accdeclnewlat(p):
  if p >= 1:
    take = 1
  else:
    compval = np.random.random(1)
    if p >= compval:
      take = 1
    else:
      take = 0
  return take
  
#The de and dm routine
def calcdedm(Lattice,Nstep,magfield,x,y,N):
  if y == 0:
    a = N-1
  else:
    a = y - 1
  if y == N-1:
    b = 0
  else:
    b = y + 1
  if x == 0:
    c = N-1
  else:
    c = x - 1
  if x == N-1:
    d = 0 
  else:
    d = x + 1
    
  sumofneighbours = Lattice[a][x] + Lattice[b][x] + Lattice[y][c] + Lattice[y][d]
  de = 2*2*Lattice[y][x]*sumofneighbours + 2*Lattice[y][x]*magfield
  dm = 2*Lattice[y][x]
  return de,dm

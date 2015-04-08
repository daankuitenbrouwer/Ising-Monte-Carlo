import numpy as np #math package
import matplotlib.pyplot as plt #Histogram package

def correlationfunction(N):#Relevant parameters

  maxdist=np.sqrt(2)*(N/2)
  distmat = np.zeros((N,N))

  nbins = int(np.sqrt(2)*(N/2)+1)
  x = np.zeros((nbins+2,1))
  binvec = np.zeros((N,1))
  dx=maxdist/nbins

#Step Zero: input matrix from other Ising file

  inputmat = np.random.uniform(-1,1,(N,N))
  inputmat = np.where(inputmat <= 0, -1, 1)

#inputmat = -1*np.ones((N,N))

  pos=[N/2,N/2]

  for j in range(N):
    for i in range(N):
      distmat[j,i]=np.sqrt((i-pos[0])**2 + (j-pos[1])**2)
      binnum=int(distmat[j,i]*nbins/(maxdist))+1
      if inputmat[j,i] == +1:
        binvec[binnum]=binvec[binnum]+1
      else:
        binvec[binnum]=binvec[binnum]-1
  
  for binnum in range(nbins+1): #Indexing to avoid divide by zero
    x[binnum+1] = (binnum+1)*dx
    binvec[binnum+1] = binvec[binnum+1]/(2*np.pi*x[binnum+1]*dx) #Divide by area of encompassing circle
    
  print binvec.shape

  return x,binvec

#Step Three: Do this for all positions in the lattice and average over them.



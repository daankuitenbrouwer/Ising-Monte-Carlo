import numpy as np #math package
import matplotlib.pyplot as plt #Histogram package
import correlationmodule
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

#Parameters:
N = 100 #size lattice
runs = N #Number of temperature points
lowT = 1 #Lower bound temperature
highT = 3 #Higher bound temperature
Trange = np.linspace(lowT,highT,runs) 
bigmat = np.zeros((N,runs))

"""
print bigmat

testmat = np.array([[[1,2],[3,4]],[[5,6],[7,8]]])

print testmat

slicemat = np.array([[9,10],[11,12]])

testmat[:,1,:] = slicemat

print testmat

"""

for run in range(runs):
  x,z = correlationmodule.correlationfunction(N)
  bigmat[:,run] = np.transpose(z)
  
print bigmat

bignum = N
X, Y = np.mgrid[:bignum, :bignum]

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
surf = ax.plot_surface(X,Y,bigmat, cmap=cm.coolwarm)
plt.show()

  
  
  
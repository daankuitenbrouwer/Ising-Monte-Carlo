import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# runs and Trange have to be specified to open the mentioned datafile

size = 22
runs = 100
Trange = np.linspace(1,3,300)

magarray22100300 = np.loadtxt("magnetizationruns%sTrange%s.txt"%(runs,len(Trange)))

size =100
runs = 30
Trange = np.linspace(1,3,40)

filee = ['magnetizationruns100Trange300.txt','magnetizationruns30Tsteps400size50.txt',
'magnetizationruns100Trange100.txt','magnetizationruns100Trange100.txt',
'magnetizationruns30Tsteps40size100.txt','magnetizationruns100Tsteps400size22lowT1highT3.txt']

# Critical temperature and exponent are defined, the theoretical plot and empirical plot are both made.
Tc = 1.859
beta = 0.5
theor_val_T = np.linspace(1.5,Tc,100)
mag = (abs(theor_val_T - Tc))**beta
plt.plot(theor_val_T,mag,color='black')
plt.plot(theor_val_T,-mag,color='black')
#plt.scatter(np.linspace(1,3,300),np.mean(np.loadtxt(filee[0]),axis=0),color='red')
#plt.scatter(np.linspace(1.2,2.6,400),np.mean(np.loadtxt(filee[1]),axis=0),color='blue')
#plt.scatter(np.linspace(1,3,100),np.mean(np.loadtxt(filee[2]),axis=0),color='pink')
#plt.scatter(np.linspace(1,3,100),np.mean(np.loadtxt(filee[3]),axis=0),color='green')
#plt.scatter(np.linspace(1,3,40),np.mean(np.loadtxt(filee[4]),axis=0),color='black')#plt.scatter(Trange,np.mean(magarray,axis=0))
plt.scatter(np.linspace(1,3,400),np.mean(np.loadtxt(filee[5]),axis=0),color='yellow')
plt.xlabel('KbT/J')
plt.ylabel('m')
#plt.text(1.8,1.3,'size = %s, runs = %s, Tsteps = %s'%(size,runs,len(Trange)))
plt.title('Critical Temperature')
#plt.text(1.8,1.2,'blue line is theoretical critical behaviour')
#plt.show()


Tlow = 1
Thigh = 3
Tsteps = 400
magarray = np.mean(np.loadtxt(filee[5]),axis=0)

pos_magarray = np.sqrt(magarray**2)
T_vals = np.linspace(Tlow,Thigh,Tsteps)
delete_lst = []
for ind,T in enumerate(T_vals):
  if T < (Tc - 0.4):
    delete_lst.append(ind)
  if T > (Tc - 0.001):
    delete_lst.append(ind)
 
delete_array = np.array(delete_lst)   
crit_T_vals = np.delete(T_vals,delete_array)
crit_pos_magarray = np.delete(pos_magarray,delete_array)

#print crit_pos_magarray,crit_T_vals
    
def beta_finder(reduced_T,beta,a):
  return a*(abs(reduced_T))**beta
 
reduced_T = crit_T_vals-Tc 
#print reduced_T
popt, pcov = curve_fit(beta_finder,reduced_T,crit_pos_magarray) #curve_fit(beta_finder,theor_val_T-Tc,mag)#
print pcov

fitted_mag = beta_finder(reduced_T,popt[0],popt[1])#beta_finder(theor_val_T-Tc,popt[0])#

plt.plot(crit_T_vals,fitted_mag,'red')#plt.plot(theor_val_T,fitted_mag,'green')#
plt.plot(np.linspace(1,3,400),pos_magarray)
plt.text(2,1.2,'beta %s'%(popt[0]))
#plt.text(2,1.2,'sigma%s'%(np.sqrt(np.diag(pcov))))
plt.text(2,1.1,'Tcrit %s'%(Tc))
plt.show()

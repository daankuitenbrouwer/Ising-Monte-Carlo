import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# runs and Trange have to be specified to open the mentioned datafile

filee = ['magnetizationruns100Trange300.txt','magnetizationruns30Tsteps400size50.txt',
'magnetizationruns100Trange100.txt','magnetizationruns100Trange100.txt',
'magnetizationruns30Tsteps40size100.txt','magnetizationruns100Tsteps400size22lowT1highT3.txt']

chosen_file = 5
Tsteps = 400

# Critical temperature and exponent are defined, the theoretical plot and empirical plot are both made.
Tc = 1.86
beta = 0.125
theor_val_T = np.linspace(1.5,Tc,Tsteps)
mag = (abs(theor_val_T - Tc))**beta
plt.plot(theor_val_T,mag,color='black')
plt.plot(theor_val_T,-mag,color='black')
plt.scatter(np.linspace(1,3,Tsteps),np.mean(abs(np.loadtxt(filee[chosen_file])),axis=0),color='yellow')
plt.xlabel('KbT/J')
plt.ylabel('m')
plt.title('Critical Temperature')


#loading and selecting the data
Tlow = 1
Thigh = 3
pos_magarray = np.mean(abs(np.loadtxt(filee[chosen_file])),axis=0)
T_vals = np.linspace(Tlow,Thigh,Tsteps)
delete_lst = []
for ind,T in enumerate(T_vals):
  if T < (Tc - 0.35):
    delete_lst.append(ind)
  if T > (Tc - 0.05):
    delete_lst.append(ind)
 
delete_array = np.array(delete_lst)   
crit_T_vals = np.delete(T_vals,delete_array)
crit_pos_magarray = np.delete(pos_magarray,delete_array)

 
#fitting for the critical exponent   
def beta_finder(reduced_T,beta,a):
  return a*(abs(reduced_T))**beta
 
reduced_T = -crit_T_vals+Tc 
popt, pcov = curve_fit(beta_finder,reduced_T,crit_pos_magarray) #curve_fit(beta_finder,theor_val_T-Tc,mag)#

delete_lst = []
for ind,T in enumerate(T_vals):
  if T < (Tc - 0.35):
    delete_lst.append(ind)
  if T > (Tc):
    delete_lst.append(ind)
 
delete_array = np.array(delete_lst)   
Plot_T_vals = np.delete(T_vals,delete_array)
reduced_plot_T = -Plot_T_vals + Tc
fitted_mag = beta_finder(reduced_plot_T,popt[0],popt[1])

plt.plot(Plot_T_vals,fitted_mag,'red')
plt.plot(np.linspace(1,3,Tsteps),pos_magarray)
plt.text(2,1.2,'beta %s'%(popt[0]))
plt.text(2,1.1,'Tcrit %s'%(Tc))
plt.show()

#plt.scatter(np.linspace(1,3,300),np.mean(np.loadtxt(filee[0]),axis=0),color='red')
#plt.scatter(np.linspace(1.2,2.6,400),np.mean(np.loadtxt(filee[1]),axis=0),color='blue')
#plt.scatter(np.linspace(1,3,100),np.mean(np.loadtxt(filee[2]),axis=0),color='pink')
#plt.scatter(np.linspace(1,3,100),np.mean(np.loadtxt(filee[3]),axis=0),color='green')
#plt.scatter(np.linspace(1,3,40),np.mean(np.loadtxt(filee[4]),axis=0),color='black')#plt.scatter(Trange,np.mean(magarray,axis=0))

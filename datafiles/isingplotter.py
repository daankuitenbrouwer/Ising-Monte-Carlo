import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# loading data
runs = 20
Tsteps = 200
size = 100
lowT = 1.3
highT = 2.7
Trange = np.linspace(lowT,highT,Tsteps)
mag_array = np.loadtxt('magnetizationruns%sTsteps%ssize%slowT%shighT%s.txt'%(runs,Tsteps,size,lowT,highT))
en_array = np.loadtxt('energyruns%sTsteps%ssize%slowT%shighT%s.txt'%(runs,Tsteps,size,lowT,highT))



# Critical temperature and exponent are defined, the theoretical plot is made.
Tc = 2.2727
beta = 0.125
theor_val_T = np.linspace(1.5,Tc,Tsteps)
theor_mag = (abs(theor_val_T - Tc))**beta
plt.plot(theor_val_T,theor_mag,color='black')
plt.plot(theor_val_T,-theor_mag,color='black')
plt.xlabel('KbT/J')
plt.ylabel('m')
plt.title('Critical Temperature')


#preparing, selecting and fitting the data for beta
pos_magarray = abs(np.mean(mag_array,axis=0))
T_vals = np.linspace(lowT,highT,Tsteps)
delete_lst = []
for ind,T in enumerate(T_vals):
  if T < (Tc - 0.3):
    delete_lst.append(ind)
  if T > (Tc -0.00001):
    delete_lst.append(ind)
 
delete_array = np.array(delete_lst)   
crit_T_vals = np.delete(T_vals,delete_array)
crit_pos_magarray = np.delete(pos_magarray,delete_array)

def beta_finder(reduced_T,beta,a):
  return a*(abs(reduced_T))**beta
 
reduced_T = -crit_T_vals+Tc 
popt, pcov = curve_fit(beta_finder,reduced_T,crit_pos_magarray) 

 
#plotting for the critical exponent beta  
delete_lst = []
for ind,T in enumerate(T_vals):
  if T < (Tc - 0.3):
    delete_lst.append(ind)
  if T > (Tc ):
    delete_lst.append(ind)
 
delete_array = np.array(delete_lst)   
Plot_T_vals = np.delete(T_vals,delete_array)
reduced_plot_T = -Plot_T_vals + Tc
fitted_mag = beta_finder(reduced_plot_T,popt[0],popt[1])

plt.plot(Plot_T_vals,fitted_mag,'red')
plt.scatter(np.linspace(lowT,highT,Tsteps),abs(np.mean(mag_array,axis=0)),color='yellow')
plt.plot(np.linspace(lowT,highT,Tsteps),pos_magarray)
plt.text(2,1.2,'beta %s'%(popt[0]))
plt.text(2,1.1,'Tcrit %s'%(Tc))
plt.show()

#preparing data for heat capacity and magnetic susceptibility
avg_m_sq_array = (np.mean(mag_array,axis=0))**2
m_sq_avg_array = np.mean(mag_array**2,axis=0)
magn_susc = abs((1./Trange)*(avg_m_sq_array - m_sq_avg_array))

avg_en_sq_array = (np.mean(en_array,axis=0))**2
en_sq_avg_array = np.mean(en_array**2,axis=0)
heat_cap = (avg_en_sq_array - en_sq_avg_array)

en_array_altern = -2*mag_array**2*size**2
avg_en_sq_array_altern = (np.mean(en_array_altern,axis=0))**2
en_sq_avg_array_altern = np.mean(en_array_altern**2,axis=0)
heat_cap_altern = (avg_en_sq_array_altern - en_sq_avg_array_altern)

plt.figure(1)
plt.subplot(121)
#plt.scatter(Trange,abs(heat_cap),color='green')
plt.scatter(Trange,abs(heat_cap_altern),color='red')
#plt.text(0,0,'heat_cap')
plt.subplot(122)
plt.scatter(Trange,magn_susc)
#plt.text(0,0,'magn_susc')
plt.title('heatcap and magnsusc')
plt.show()



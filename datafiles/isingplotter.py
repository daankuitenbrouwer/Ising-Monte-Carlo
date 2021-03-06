import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# loading data
runs = 20
Tsteps = 200
size = 100
lowT = 1.3
highT = 2.7
Tc = 2.2727
Trange = np.linspace(lowT,highT,Tsteps)
mag_array = np.loadtxt('magnetizationruns%sTsteps%ssize%slowT%shighT%s.txt'%(runs,Tsteps,size,lowT,highT))
#en_array = np.loadtxt('energyruns%sTsteps%ssize%slowT%shighT%s.txt'%(runs,Tsteps,size,lowT,highT))


def beta(Tc,Tsteps,size,lowT,highT,Trange,mag_array,en_array):
    # Critical temperature and exponent are defined, the theoretical plot is made.
    beta = 0.125
    theor_val_T = np.linspace(1.5,Tc,Tsteps)
    theor_mag = (abs(theor_val_T - Tc))**beta
    #plt.plot(theor_val_T,theor_mag,color='black')
    #plt.plot(theor_val_T,-theor_mag,color='black')
    plt.xlabel('KbT/J')
    plt.ylabel('m')
    plt.title('Critical Temperature')


    #preparing, selecting and fitting the data for beta
    pos_magarray = abs(np.mean(mag_array,axis=0))
    Trange = np.linspace(lowT,highT,Tsteps)
    delete_lst = []
    for ind,T in enumerate(Trange):
      if T < (Tc - 0.3):
        delete_lst.append(ind)
      if T > (Tc -0.00001):
        delete_lst.append(ind)
     
    crit_Trange = np.delete(Trange,np.array(delete_lst))
    crit_pos_magarray = np.delete(pos_magarray,np.array(delete_lst))

    def beta_finder(reduced_T,beta,a):
      return a*(abs(reduced_T))**beta
     
    reduced_T = -crit_Trange+Tc 
    popt, pcov = curve_fit(beta_finder,reduced_T,crit_pos_magarray) 

     
    #plotting for the critical exponent beta  
    delete_lst = []
    for ind,T in enumerate(Trange):
      if T < (Tc - 0.3):
        delete_lst.append(ind)
      if T > (Tc ):
        delete_lst.append(ind)
     
    Plot_Trange = np.delete(Trange,np.array(delete_lst))
    reduced_plot_T = -Plot_Trange + Tc
    fitted_mag = beta_finder(reduced_plot_T,popt[0],popt[1])
    
    print pcov

    plt.plot(Plot_Trange,fitted_mag,'red')
    plt.scatter(np.linspace(lowT,highT,Tsteps),abs(np.mean(mag_array,axis=0)),color='green')
    #plt.plot(np.linspace(lowT,highT,Tsteps),pos_magarray)
    plt.text(2.4,1.1,r'$\beta$ = %s'%(round(popt[0],4)))
    plt.text(2.4,1.05,'Tcrit = %s'%(Tc))
    plt.show()
    return

def magn_susceptibility(Tc,Tsteps,size,lowT,highT,Trange,mag_array):
    #preparing data for heat capacity 
    avg_m_sq_array = (np.mean(mag_array,axis=0))**2
    m_sq_avg_array = np.mean(mag_array**2,axis=0)
    magn_susc = abs((1./Trange)*(avg_m_sq_array - m_sq_avg_array))

    #fitting for chi_finder
    delete_lst_chi = []
    for ind,T in enumerate(Trange):
      if T < (Tc + 0.02 ):
        delete_lst_chi.append(ind)
      if T > (Tc + 0.4):
        delete_lst_chi.append(ind)
    crit_Trange_chi = np.delete(Trange,np.array(delete_lst_chi))
    crit_pos_magarray_chi = np.delete(magn_susc,np.array(delete_lst_chi))

    #print crit_Trange_chi,crit_pos_magarray_chi

    def gamma_finder(reduced_T,gamma,a):
      return a*reduced_T**-gamma
      
    reduced_T_chi = abs(crit_Trange_chi-Tc)
    popt_gamma, pcov_gamma = curve_fit(gamma_finder,reduced_T_chi,crit_pos_magarray_chi)

    print popt_gamma,pcov_gamma

    #Plotting chi
    fitted_mag_chi = gamma_finder(reduced_T_chi,popt_gamma[0],popt_gamma[1])
    theor_mag_chi = gamma_finder(abs(Trange - Tc),(7./4.),1./500000)
    plt.plot(crit_Trange_chi,fitted_mag_chi,'red')
    #plt.plot(Trange,theor_mag_chi,'black')
    plt.scatter(Trange,magn_susc,color='green')
    #plt.plot(np.linspace(lowT,highT,Tsteps),magn_susc)
    plt.text(1.3,0.105,r'$\gamma$ = %s'%(round(popt_gamma[0],4)))
    plt.text(1.3,0.1,'Tcrit %s'%(Tc))
    plt.xlabel('KbT/J')
    plt.ylabel(r'$\chi$')
    plt.title('Magnetic Susceptibility')
    plt.show()
    return

def heat_capacity(Tc,Tsteps,size,lowT,highT,Trange,mag_array,en_array):
    #preparing data for heat capacity
    avg_en_sq_array = (np.mean(en_array,axis=0))**2
    en_sq_avg_array = np.mean(en_array**2,axis=0)
    heat_cap = abs((1./Trange**2)*(avg_en_sq_array - en_sq_avg_array))

    en_array_altern = -2*mag_array**2*size**2
    avg_en_sq_array_altern = (np.mean(en_array_altern,axis=0))**2
    en_sq_avg_array_altern = np.mean(en_array_altern**2,axis=0)
    heat_cap_altern = abs((1./Trange**2)*(avg_en_sq_array_altern - en_sq_avg_array_altern))


    #fitting for Ch
    delete_lst_Ch = []
    for ind,T in enumerate(Trange):
      if T < (Tc - 0.2):
        delete_lst_Ch.append(ind)
      if T > (Tc -0.00001):
        delete_lst_Ch.append(ind)
    crit_Trange_Ch = np.delete(Trange,np.array(delete_lst_Ch))
    crit_pos_magarray_Ch = np.delete(heat_cap_altern,np.array(delete_lst_Ch))

    '''

    def alfa_finder(reduced_T,alfa,a):
      return a*reduced_T**-alfa
      
    reduced_T_Ch = abs(crit_Trange_Ch-Tc)
    popt_alfa, pcov_alfa = curve_fit(alfa_finder,reduced_T_Ch,abs(crit_pos_magarray_Ch))
     '''
     
    def Ch_fit(crit_Trange,a,b,B):
      return a + b*crit_Trange_Ch**B
      
    #popt_Ch,pcov_Ch = curve_fit(Ch_fit,crit_Trange_Ch,abs(crit_pos_magarray_Ch))

    #Plotting Ch
    #fitted_mag_Ch = Ch_fit(crit_Trange_Ch,popt_Ch[0],popt_Ch[1],popt_Ch[2])
    #theor_mag_Ch = alfa_finder(reduced_T_Ch,0.1,100000)
    #plt.plot(crit_Trange_Ch,theor_mag_Ch,'green')
    #plt.plot(crit_Trange_Ch,fitted_mag_Ch,'red')
    plt.scatter(np.linspace(lowT,highT,Tsteps),heat_cap_altern,color='green')
    #plt.plot(np.linspace(lowT,highT,Tsteps),heat_cap)
    #plt.text(2,100002,'a=%s, b=%s, B=%s'%(popt_Ch[0],popt_Ch[1],popt_Ch[2]))
    #plt.text(2,1.1,'Tcrit %s'%(Tc))
    plt.title('Heat Capacity')
    plt.xlabel('KbT/J')
    plt.ylabel('Ch')
    plt.show()
    return

def mag_var_h():
    # loading data
    runs = 7
    runs2 = 10
    runs3 = 10
    hsteps = 600
    hsteps2 = 300
    hsteps3 = 100
    hsteps4 = 100
    size = 100
    lowh = 0
    highh = 0.12
    highh2 = 0.2
    highh3 = 0.2
    highh4 = 0.2
    hrange = np.linspace(lowh,highh,hsteps)
    hrange2 = np.linspace(lowh,highh2,hsteps2)
    hrange3 = np.linspace(lowh,highh3,hsteps3)
    hrange4 = np.linspace(lowh,highh3,hsteps3)
    mag_h_array = np.mean(abs(np.loadtxt('magnetizationruns%shsteps%ssize%slowh%shighh%s.txt'%(runs,hsteps,size,lowh,highh))),axis=0)
    mag_h_array2 = np.mean(abs(np.loadtxt('magnetizationruns%shsteps%ssize%slowh%shighh%s.txt'%(runs2,hsteps2,size,lowh,highh2))),axis=0)
    mag_h_array3 = np.mean(abs(np.loadtxt('magnetizationruns%shsteps%ssize%slowh%shighh%s.txt'%(runs3,hsteps3,size,lowh,highh3))),axis=0)
    mag_h_array4 = np.mean(abs(np.loadtxt('magnetizationruns%shsteps%ssize%slowh%shighh%sNiter1000.txt'%(runs3,hsteps3,size,lowh,highh3))),axis=0)

    #preparing data for delta finder
    delete_lst_delta = []
    for ind,h in enumerate(hrange):
      if h < 0.10:
        delete_lst_delta.append(ind)
      if h > 0.14:
        delete_lst_delta.append(ind)
    crit_h_vals_delta = np.delete(hrange,np.array(delete_lst_delta))
    crit_mag_h_array = np.delete(mag_h_array,np.array(delete_lst_delta))

    def delta_finder(h,b,invdelta):
      return  b*h**(invdelta)
      
    popt_delta, pcov_delta = curve_fit(delta_finder,crit_h_vals_delta,crit_mag_h_array)
    fitted_mag_delta = delta_finder(crit_h_vals_delta,popt_delta[0],popt_delta[1])
    
    theor_h_vals = np.linspace(0,0.2,200)
    theor_mag_delta = delta_finder(theor_h_vals,1.,1./15.)

    #plt.plot(crit_h_vals_delta,fitted_mag_delta,'red')
    plt.plot(theor_h_vals,theor_mag_delta,'blue')
    plt.scatter(hrange,mag_h_array,color='green')
    plt.scatter(hrange2,mag_h_array2,color='grey')
    plt.scatter(hrange3,mag_h_array3,color='black')
    plt.scatter(hrange4,mag_h_array4,color='yellow')
    #plt.text(0,1,'delta%s'%(1./popt_delta[1]))
    plt.xlabel('h')
    plt.ylabel('m')
    plt.title('variation magnetization for small magnetic field at Tc')
    plt.show()
    return


#heat_capacity(Tc,Tsteps,size,lowT,highT,Trange,mag_array,en_array)
mag_var_h()
#beta(Tc,Tsteps,size,lowT,highT,Trange,mag_array,en_array)
#magn_susceptibility(Tc,Tsteps,size,lowT,highT,Trange,mag_array)

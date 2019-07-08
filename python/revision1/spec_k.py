import matplotlib.pyplot as plt
import numpy as np
import os

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

m_list = ['1','4','8'] 
eqn = 'BO'
nk = 21
ip = 50

Bu_max = 0.25



fig, ax = plt.subplots(1)

for m in m_list:

    path_spec = scratch_location+folder+eqn+'/m'+m+'_U0.25_n256/output/h9w.dat'
            
    if os.path.isfile(path_spec):

        spec = np.loadtxt(path_spec)      #spectrum at all times
        spec = spec[nk*ip:nk*(ip+1),:]    #only take the spectrum at the time of interest, t = ip inertial periods
        tot_energy = np.sum(spec[:,1])    #Compute the total energy at that level
        spec[:,1] = spec[:,1]/tot_energy  #Normalize 

        Bu = Bu_max/(int(m)**2)  #to be multiplied by k^2
#        plt.plot(Bu*spec[:,0]**2,spec[:,1],"-*",ms=8,label="BQ , m'="+m)
        plt.plot(spec[:,0],spec[:,1],"-*",ms=8,label="BQ , m'="+m)

ax.grid(color='k', linestyle='-', linewidth=0.1)
plt.yscale("log")
plt.xscale("log")
#plt.title('Energy distribution after 50 inertial periods')
#plt.legend(loc='lower left')        
plt.legend(loc='upper right')        
#plt.xlabel("Bu(k',m')")
plt.xlabel(r"$k_h$")
plt.ylabel('Normalized Energy Density',rotation=90,labelpad=10)
plt.xlim(right=30)
#plt.xlim(right=100)
#plt.ylim((1e-5,1e-0))
plt.axvline(x=1,linewidth=0.5, color='k')
plt.show()
#plt.savefig('plots/spec_k.eps')


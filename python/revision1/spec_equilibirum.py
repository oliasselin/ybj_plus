import matplotlib.pyplot as plt
import numpy as np
import os

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/revision1/'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

m = '8' 
eqn = 'YBJp'
nk = 85#21 #Number of kh's in the spectrum
ip_min = 280
ip_max = 320

diss = '_c0.01'

fig, ax = plt.subplots(1)

path_spec        = scratch_location+folder+eqn+'/m'+m+'_U0.25_n256_hres'+diss+'/output/h0w.dat'
path_spec_nodisp = scratch_location+folder+eqn+'/m'+m+'_U0.25_n256_hres_nodisp'+diss+'/output/h0w.dat'
            
if os.path.isfile(path_spec):

    spec        = np.loadtxt(path_spec)      #spectrum at all times
    spec_nodisp = np.loadtxt(path_spec_nodisp)      #spectrum at all times

    spec_ip = np.zeros((nk,4,ip_max-ip_min))
    spec_nodisp_ip = np.zeros((nk,4,ip_max-ip_min))

    for ip in range(ip_min,ip_max):

        ipp=ip-ip_min
        spec_ip[:,:,ipp]        = spec[nk*ip:nk*(ip+1),:]    #only take the spectrum at the time of interest, t = ip inertial periods
        spec_nodisp_ip[:,:,ipp] = spec_nodisp[nk*ip:nk*(ip+1),:]    #only take the spectrum at the time of interest, t = ip inertial periods



    ave_spec        = np.mean(spec_ip,2)
    ave_spec_nodisp = np.mean(spec_nodisp_ip,2)

    tot_energy = np.sum(spec_ip[:,1,0])

    plt.plot(ave_spec[:,0],ave_spec[:,1]/tot_energy,"-",color='mediumorchid',lw=2.,ms=8,label="YBJ$^+$, m'="+m+", with dispersion")
    plt.plot(ave_spec_nodisp[:,0],ave_spec_nodisp[:,1]/tot_energy,"-",color='olive',lw=2.,ms=8,label="YBJ$^+$, m'="+m+", no dispersion")

ax.grid(color='k', linestyle='-', linewidth=0.1)
plt.yscale("log")
plt.xscale("log")
#plt.title('Energy distribution after 50 inertial periods')
plt.legend(loc='best',fontsize='small')        
#plt.legend(loc='upper right')        
#plt.xlabel("Bu(k',m')")
plt.xlabel(r"$k_h$")
plt.ylabel('Normalized Energy Density',rotation=90,labelpad=10)
#plt.xlim(right=10)
#plt.xlim(right=2)
plt.ylim((1.e-6,1))
#plt.axvline(x=1,linewidth=0.5, color='k')
#plt.show()
plt.savefig('plots/spec_eq.eps')


import matplotlib.pyplot as plt
import numpy as np
import os

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

m_list = ['1']#['1','4','8'] 
eqn = 'YBJp'
nk = 21
ip_max = 51
ip_list=[1,10,20,35,50]

Bu_max = 0.25


vres = 256
dz   = 2.*np.pi/vres
U_scale = 0.25
Uw_scale = 2.5e-5

for m in m_list:

    path_spec = scratch_location+folder+eqn+'/m'+m+'_U0.25_n256/output/h0w.dat'

    pred_energy= 0.5*((Uw_scale/U_scale)**2)*(np.cos(int(m)*(2*np.pi-dz/2.)/2.))**2
            
    if os.path.isfile(path_spec):

        spec = np.loadtxt(path_spec)      #spectrum at all times

        time = np.zeros(ip_max)
        tot_energy = np.zeros(ip_max)
        cumulative_energy = np.zeros((ip_max,nk))

        for ip in range(ip_max):

            time[ip]=ip
            spec_ip = spec[nk*ip:nk*(ip+1),:]    #only take the spectrum at the time of interest, t = ip inertial periods

            #Modification of the output: the kh > 0 modes should be multiplied by 2 (not the kh=0 mode)
            spec_ip[1:,1] = spec_ip[1:,1]*2.
            #############################

            tot_energy[ip] = np.sum(spec_ip[:,1])    #Compute the total energy at that level

            
            for cumul in range(nk):
                cumulative_energy[ip,cumul] = np.sum(spec_ip[0:cumul+1,1])
         
            Bu = Bu_max*(spec_ip[:,0]**2)/(int(m)**2)  

        for ip in ip_list:    
            plt.plot(Bu,cumulative_energy[ip,:]/pred_energy,"-*",ms=8,label="YBJ$^+$, ip="+str(ip)+", m'="+m)

    plt.xlim(right=1)
    plt.legend(loc='lower right')
#    plt.plot(time,tot_energy/pred_energy)#
    plt.xlabel(r"$Bu = (Nk/fm)^2$")
    plt.ylabel('Cumulative Energy',rotation=90,labelpad=10)
    plt.show()
#        print "Integrated spectrum=",tot_energy,"Prediction",pred_energy


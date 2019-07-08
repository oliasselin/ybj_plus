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
ip_max = 50

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

        for ip in range(ip_max):

            time[ip]=ip
            spec_ip = spec[nk*ip:nk*(ip+1),:]    #only take the spectrum at the time of interest, t = ip inertial periods
            tot_energy[ip] = np.sum(spec_ip[:,1])    #Compute the total energy at that level
            #spec[:,1] = spec[:,1]/tot_energy  #Normalize 



    plt.plot(time,tot_energy/pred_energy)
    plt.show()
#        print "Integrated spectrum=",tot_energy,"Prediction",pred_energy


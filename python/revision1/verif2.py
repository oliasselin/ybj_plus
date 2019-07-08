import matplotlib.pyplot as plt
import numpy as np
import os

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

m_list = ['8']#['1','4','8'] 
eqn = 'YBJ'
nk = 21
ip_max = 50

Bu_max = 0.25


vres = 256
dz   = 2.*np.pi/vres
U_scale = 0.25
Uw_scale = 2.5e-5
u = str(U_scale)


for m in m_list:

    time = np.zeros(ip_max)
    tot_energy = np.zeros(ip_max)
    
    for ts in range(ip_max):

        run = 'm'+m+'_U'+u+'_n'+str(vres)
        spaces_ts = (3-len(str(ts)))*' '

        path_ybj  = scratch_location+folder+eqn+'/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'

        pred_energy= 0.5*((Uw_scale/U_scale)**2)*(np.cos(int(m)*(2*np.pi-dz/2.)/2.))**2
            
        if os.path.isfile(path_ybj):

            wke_slice = np.loadtxt(path_ybj)      #spectrum at all times
            tot_energy[ts]= np.average(wke_slice)
            time[ts] = ts

    plt.plot(time,tot_energy/pred_energy)
    plt.show()
#        print "Integrated spectrum=",tot_energy,"Prediction",pred_energy


import matplotlib.pyplot as plt
import numpy as np
import os

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

m_list = ['1','4','8','12'] 
nm = len(m_list)
eqn = 'YBJp'
nk = 21
ip = 50

ymin=80
ymax=100

#Plot vertical lines at various frequencies?
plot_vlines=1
freq_list = ['1.1','1.2','1.3','1.4']

Bu_max = 0.25

vres = 256
dz   = 2.*np.pi/vres
U_scale = 0.25
Uw_scale = 2.5e-5

fig = plt.figure()
ax = fig.add_subplot(111)

for im,m in enumerate(m_list):

    path_spec = scratch_location+folder+eqn+'/m'+m+'_U0.25_n256/output/h0w.dat'
    pred_energy= 0.5*((Uw_scale/U_scale)**2)*(np.cos(int(m)*(2*np.pi-dz/2.)/2.))**2
            
    if os.path.isfile(path_spec):

        spec = np.loadtxt(path_spec)      #spectrum at all times
        cumulative_energy = np.zeros((nm,nk))

        spec_ip = spec[nk*ip:nk*(ip+1),:]    #only take the spectrum at the time of interest, t = ip inertial periods
        #Modification of the output: the kh > 0 modes should be multiplied by 2 (not the kh=0 mode)
        spec_ip[1:,1] = spec_ip[1:,1]*2.
        ################################
            
        for cumul in range(nk):
            cumulative_energy[im,cumul] = np.sum(spec_ip[0:cumul+1,1])
         
        Bu = Bu_max*(spec_ip[:,0]**2)/(int(m)**2)  

#        plt.plot(Bu,100*cumulative_energy[im,:]/pred_energy,"-*",ms=8,label="YBJ$^+$, ip="+str(ip)+", m'="+m)
#        plt.plot(Bu,100*cumulative_energy[im,:]/pred_energy,"-+",linewidth=0.25,ms=10,label="YBJ$^+$, ip="+str(ip)+", m'="+m)
        plt.plot(Bu,100*cumulative_energy[im,:]/pred_energy,"-*",linewidth=1.,label="YBJ$^+$, $Ro$ = 0.05, $m'$="+m)

plt.xlim(right=1)
plt.ylim((ymin,ymax))
plt.legend(loc='lower right')
plt.xlabel(r"$Bu = (Nk/fm)^2$")
plt.ylabel('Cumulative Energy (%)',rotation=90,labelpad=10)
#plt.xaxis.grid(True, which='minor')
#plt.minorticks_on()
ax.yaxis.grid(True, which='minor')




plt.grid(color='k', linestyle='-', linewidth=0.05)


if(plot_vlines==1):
    #Plot vertical lines at various frequencies
    for freq in freq_list:
        
#        ff = float(freq)-1
#        Bu_coord = 4*(float(freq)-1))/(2-(float(freq)-1))
        Bu_coord = 4*(float(freq)-1)/(2-(float(freq)-1))

        plt.axvline(x=Bu_coord, linewidth=1., color='k')
        plt.text(Bu_coord-0.05, 101,r'$\omega='+freq+'f$')

#plt.show()
plt.savefig('plots/cumul_m.eps')


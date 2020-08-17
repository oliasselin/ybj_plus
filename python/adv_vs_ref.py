import matplotlib.pyplot as plt
import numpy as np
import os

show=0

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'gamma/dipole_2xU'#'N2_1e5_dla'

location = scratch_location+folder+run
gamma = np.loadtxt(scratch_location+folder+run+'/output/gamma.dat')

if not os.path.exists('plots/'+run+'/'):
    os.makedirs('plots/'+run+'/')

plt.plot(gamma[:,0],gamma[:,1],'-g',label=r'$\Gamma_a$',linewidth=2.)
plt.plot(gamma[:,0],gamma[:,2],'-b',label=r"$\Gamma_r$",linewidth=2.)
plt.plot(gamma[:,0],gamma[:,3],'-k',label=r"$\Gamma_d$",linewidth=2.)
plt.plot(gamma[:,0],gamma[:,0]*0,'--k',label='_nolegend_',linewidth=0.5)


#plt.ylim((0.,0.0015))
plt.title(r"WPE production, run ="+run)

plt.legend(loc='upper left',prop={'size': 9})
plt.xlabel(r'$t$ (days)')
plt.ylabel(r'WPE production (W/kg)')
    
if(show==1):
    plt.show()
else:
    plt.savefig('plots/'+run+'/gamma.eps',bbox_inches='tight')


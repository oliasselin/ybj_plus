import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os
import subprocess
import sys

timestep=0.1 #0.1 #Fraction of an inertial period between slices

#ts_plot=[30]

hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'N2_1e5'

plot_slice=1
colormap='RdBu_r' 
ncases=1
nts=1

aspect=0.6
depth = 1000  #m
Dz= 3000
gp_depth = int(vres*depth/Dz)
dz= Dz/vres

nyticks=4
tylabels=np.arange(0,depth+1,depth/nyticks)
ticksy_loc=np.arange(0,gp_depth,(gp_depth-1.)/nyticks)

xyrange_km = 35#np.sqrt(2)*25
xrange_km  = xyrange_km*np.cos(np.deg2rad(45))
Dx_km=100
gp_del= int(round(hres*xrange_km/Dx_km))
x0 = int(hres/2-gp_del)
x1 = int(hres/2+gp_del)
gp_xrange=x1-x0


nxticks=2
#txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)
###


omega_min=-0.3
omega_max=0.3
wke_threshold = 1e-6
wke_levels=[1e-6, 1e-5, 1e-4, 1e-3, 1e-2]

#Time range (time should be in fractions of inertial period. Fraction = timestep var above)                                                                                                                   
ts_min=0
ts_max=203


#Create folder for plots if it doesn't exists                                                                                                                                            
if not os.path.exists('plots/'+run+'/omega/'):
    os.makedirs('plots/'+run+'/omega/')
#Suboptimal: load vertical slice                                                                                                                              


#t=0: load three times
ts=0

spaces_ts = (3-len(str(ts)))*' '
path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'
         

#Load the data file                                                                                                                                                                  
if os.path.isfile(path_lar):

    f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                                             
    f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                                                
    
    g_lar   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                            
    g_lai   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                       
    
    lar_m = g_lar[0:gp_depth,x0:x1]
    lai_m = g_lai[0:gp_depth,x0:x1]
    
ts=1

spaces_ts = (3-len(str(ts)))*' '
path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'

#Load the data file                                                                                                                                                                  
if os.path.isfile(path_lar):

    f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                                             
    f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                                                
    
    g_lar   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                            
    g_lai   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                       
    
    lar_0 = g_lar[0:gp_depth,x0:x1]
    lai_0 = g_lai[0:gp_depth,x0:x1]

    
for ts in range(ts_min+2,ts_max):

    
    ts0 = ts-1
    print ts0

    spaces_ts = (3-len(str(ts)))*' '
    path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
    path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'


    #Load the data file                                                                                                                                                      
    if os.path.isfile(path_lar):

        time=1.*ts0*timestep #time in inertial periods         

        f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                                             
        f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                                                   

        g_lar   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                   
        g_lai   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                         

        lar_p = g_lar[0:gp_depth,x0:x1]                                                                                                                                
        lai_p = g_lai[0:gp_depth,x0:x1]

        #Formula is sigma = 1/2dt (blah blah). If dt is given in inertial periods and we want sig/f, then the formula is 1/4pi dt' (blah blah)

        omega = -(1./(4.*np.pi*timestep))*(  (lai_p-lai_m)*lar_0 - (lar_p-lar_m)*lai_0 )/(lar_0*lar_0 + lai_0*lai_0)

        nan_indices = np.where(lar_0*lar_0 + lai_0*lai_0 < 2.*wke_threshold)
        omega[nan_indices] = 0.

        zeros_ts = (3-len(str(ts0)))*'0'
        np.savetxt('plots/'+run+'/omega/omega_t'+zeros_ts+str(ts0)+'.dat',omega)


        if(plot_slice==1):# and (ts in ts_plot)):

            print "Making a plot of omega..."

            #Produce the ncases x nts multiplot
            fig = plt.figure(figsize=(6, 4))                        
            grid = AxesGrid(fig, 111,
                            nrows_ncols=(ncases,nts),
                            axes_pad=0.05,
                            cbar_mode='single',
                            cbar_location='right',
                            cbar_pad=0.1
                            )
    
            for idx,ax in enumerate(grid):
        
                ax.get_xaxis().set_ticks(ticksx_loc)
                ax.set_xticklabels(txlabels)
                ax.get_xaxis().set_label('Depth (m)')
                
                ax.get_yaxis().set_ticks(ticksy_loc)
                ax.set_yticklabels(tylabels)
                ax.get_yaxis().set_label('Depth (m)')

                im = ax.imshow(omega,cmap=colormap,aspect=aspect,vmin=omega_min,vmax=omega_max)
                ax.contour(0.5*(lar_0*lar_0+lai_0*lai_0),levels=wke_levels,colors='k')#, levels, colors='k', origin='upper', extent=extent)
        
                cbar = ax.cax.colorbar(im)# ticks=[1e-8, 1e-2])
                cbar = grid.cbar_axes[0].colorbar(im)

                time_title = '%.1f' % time        
                ax.set_title(r'$\sigma/f$, $t =$ '+time_title+' inertial periods',fontsize=12)
                ax.text(-15, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
                ax.text(gp_del, gp_depth+20,r"$x'$ (km)",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)

#                ax.plot(np.arange(1,2*gp_del), np.arange(1,2*gp_del), color='y')

#                plt.show()
                zeros_ts = (3-len(str(ts0)))*'0'
                plt.savefig('plots/'+run+'/omega/omega_t'+zeros_ts+str(ts0)+'.png',bbox_inches='tight')
    
        lar_m = lar_0
        lai_m = lai_0

        lar_0 = lar_p
        lai_0 = lai_p



make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/omega/*.png plots/'+run+'/omega/omega.gif'
p = subprocess.Popen(make_gif, shell = True)
os.waitpid(p.pid, 0)

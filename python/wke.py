import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os
import subprocess
import sys

timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'attempt3_ro50'

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


wke_min=1e-8
wke_max=1e-2
ticks=np.logspace(-8,-2,7)


#Time range (time should be in fractions of inertial period. Fraction = timestep var above)                                                                                                                   
ts_min=0
ts_max=201


#Create folder for plots if it doesn't exists                                                                                                                                            
if not os.path.exists('plots/'+run+'/wke/'):
    os.makedirs('plots/'+run+'/wke/')
#Suboptimal: load vertical slice                                                                                                                              

    
for ts in range(ts_min,ts_max):

    
    print "ts=",ts

    spaces_ts = (3-len(str(ts)))*' '
    path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
    path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'


    #Load the data file                                                                                                                                                      
    if os.path.isfile(path_lar):

        time=ts*timestep #time in inertial periods         

        f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                                             
        f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                                                   

        g_lar   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                   
        g_lai   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                         

        lar = g_lar[0:gp_depth,x0:x1]                                                                                                                                
        lai = g_lai[0:gp_depth,x0:x1]

        wke = 0.5*(lar*lar+lai*lai)

        if(plot_slice==1):

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

                im = ax.imshow(wke,cmap=colormap,norm=colors.LogNorm(vmin=wke_min, vmax=wke_max),aspect=aspect)

                cbar = ax.cax.colorbar(im, ticks=ticks)
                cbar.ax.set_yticklabels(ticks)
                


                #im = ax.imshow(omega,cmap=colormap,aspect=aspect,vmin=omega_min,vmax=omega_max)
                #ax.contour(0.5*(lar_0*lar_0+lai_0*lai_0),levels=wke_levels,colors='k')#, levels, colors='k', origin='upper', extent=extent)
        
                #cbar = ax.cax.colorbar(im)# ticks=[1e-8, 1e-2])
                #cbar = grid.cbar_axes[0].colorbar(im)

                time_title = '%.1f' % time        
                ax.set_title(r'WKE (m/s)$^2$, $t =$ '+time_title+' inertial periods',fontsize=12)
                #ax.set_title(r'$\omega/f$, $t =$ '+time_title+' inertial periods',fontsize=12)
                ax.text(-15, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
                ax.text(gp_del, gp_depth+20,r'$x_{cs}$ (km)',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)


#                plt.show()
                zeros_ts = (3-len(str(ts)))*'0'
                plt.savefig('plots/'+run+'/wke/wke_t'+zeros_ts+str(ts)+'.png',bbox_inches='tight')


make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/wke/*.png plots/'+run+'/wke/wke.gif'
p = subprocess.Popen(make_gif, shell = True)
os.waitpid(p.pid, 0)

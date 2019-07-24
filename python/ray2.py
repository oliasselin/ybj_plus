import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os
import subprocess
import sys

timestep=0.1 #0.1 #Fraction of an inertial period between slices
ts_plot=np.arange(0,201,1)#[50,100,150,200]





hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'attempt3_ro50'
N2=1e-5
f2=np.power(1.24e-4,2)
zeta_max=2*np.pi/(100000*1.24e-4)

print(zeta_max)

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
Dx=Dx_km*1000
gp_del= int(round(hres*xrange_km/Dx_km))
x0 = int(hres/2-gp_del)
x1 = int(hres/2+gp_del)
gp_xrange=x1-x0


xstart_list=(int(1*gp_del/4),int(gp_del/2),int(3*gp_del/4),gp_del,int(5*gp_del/4),int(3*gp_del/2),int(7*gp_del/4),)
dx_cs = 1000.*xyrange_km/gp_del   #dx_cs in meters

nxticks=2
#txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)
###


omega_min=-0.3
omega_max=0.3
wke_threshold = 1e-6
wke_levels=[1e-6, 1e-5, 1e-4, 1e-3, 1e-2]


#Create folder for plots if it doesn't exists                                                                                                                 
if not os.path.exists('plots/'+run+'/ray/'):
    os.makedirs('plots/'+run+'/ray/')

for ts in ts_plot:

    print "ts=",ts

    zeros_ts = (3-len(str(ts)))*'0'
    path_omega  = 'plots/'+run+'/omega/omega_t'+zeros_ts+str(ts)+'.dat'         

    #Load the data file                                                                                                                                                      
    if os.path.isfile(path_omega):

        omega = np.loadtxt(path_omega)

        z_traj  = np.zeros((len(xstart_list),omega.shape[1]))
        iz_traj = np.zeros((len(xstart_list),omega.shape[1]))
        ix_traj = np.zeros((len(xstart_list),omega.shape[1]))


        #Compute trajectories
        for traj,xstart in enumerate(xstart_list):
         
            iz_traj[traj,0]=1
            ix_traj[traj,0]=xstart

            z_traj[traj,0]=-iz_traj[traj,0]*dz

            for it in range(0,omega.shape[1]-xstart-1):

                ix_traj[traj,it] = xstart + it
                x_cs = -xyrange_km*1000+ix_traj[traj,it]*dx_cs   #x_cs in m

                #Compute slope
                zeta_eff=-zeta_max*np.sin(2.*np.pi*x_cs/(np.sqrt(2)*Dx))
                feff2=np.power( 1+0.5*zeta_eff, 2)
                sig  = 1.+omega[iz_traj[traj,it],ix_traj[traj,it]]
                sig2 = np.power(sig,2)

#            print "sig2,feff,N2",sig2,feff2,N2

                arg= f2*(sig2-feff2)/N2
                if(arg>0):
                    slope=-np.sqrt(f2*(sig2-feff2)/N2)
                else:
                    slope=0.
                    print "feff2>sig2: break loop"
#                    ix_traj[traj,it:]=ix_traj[traj,it]
#                    iz_traj[traj,it:]=iz_traj[traj,it]
                    break

                


                z_traj[traj,it+1] = z_traj[traj,it] + slope*dx_cs

                print "slope,z,x=",slope,z_traj[traj,it+1],x_cs

                #Transpose to grid points
                iz_traj[traj,it+1]=round(-z_traj[traj,it+1]/dz)
                if(iz_traj[traj,it+1]>=gp_depth):
                    iz_traj[traj,it+1]=gp_depth-1
                    print "reached the bottom at x_cs=",x_cs
                
            ix_traj[traj,it:]=ix_traj[traj,it]
            iz_traj[traj,it:]=iz_traj[traj,it]        
        

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
                
                im = ax.imshow(omega,cmap=colormap,aspect=aspect,vmin=omega_min,vmax=omega_max)
                cbar = ax.cax.colorbar(im)# ticks=[1e-8, 1e-2])                                                                                                          
                cbar = grid.cbar_axes[0].colorbar(im)
                
                time=ts*timestep
                time_title = '%.1f' % time
                ax.set_title(r'$\omega/f$, $t =$ '+time_title+' inertial periods',fontsize=12)
                ax.text(-15, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
                ax.text(gp_del, gp_depth+20,r'$x_{cs}$ (km)',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)
                
                for traj,xstart in enumerate(xstart_list):
                    ax.plot(ix_traj[traj,:],iz_traj[traj,:],'-k',linewidth=2.5)                                                                    
                
#                plt.show()                                                                                                                                                              
            zeros_ts = (3-len(str(ts)))*'0'
            plt.savefig('plots/'+run+'/ray/ray_t'+zeros_ts+str(ts)+'.png',bbox_inches='tight')

make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/ray/*.png plots/'+run+'/ray/ray.gif'
p = subprocess.Popen(make_gif, shell = True)
os.waitpid(p.pid, 0)

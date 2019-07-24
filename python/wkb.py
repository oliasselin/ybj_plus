import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os
import subprocess
import sys

ts_plot=[5] #np.arange(1,202,1)#[50,100,150,200]

field = 6

field_title = ['$\sigma/(\zeta_{max}/2)$','$Bu/\zeta_{max}$','Doppler shift $/(\zeta_{max}/2)$','$R_{WKB}$','$(\sigma-f_{eff})/(\zeta_{max}/2)$','$f_{eff}$']
field_name  = ['sig','bu','ds','res','smz','feff']









#Run specification##############################################
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'N2_1e5'

N2=1e-5
cor = 1.24e-4
f2=np.power(cor,2)
zeta_max=2*np.pi/(100000*1.24e-4)  #Simply = k/f
U_scale = 0.5

hres=256
vres=256
timestep=0.1 #0.1 #Fraction of an inertial period between slices
################################################################



depth = 1000  #m
Dz= 3000
gp_depth = int(vres*depth/Dz)
dz= Dz/vres

xyrange_km = 35#np.sqrt(2)*25
xrange_km  = xyrange_km*np.cos(np.deg2rad(45))
Dx_km=100
Dx=Dx_km*1000
gp_del= int(round(hres*xrange_km/Dx_km))
x0 = int(hres/2-gp_del)
x1 = int(hres/2+gp_del)
gp_xrange=x1-x0
dx = xyrange_km*1000.*2./(1.*gp_xrange)

gpx = np.arange(0,gp_xrange,1)
gpz = np.arange(0,gp_depth,1)

xp_m = -xyrange_km*1000 + dx*gpx

k_over_sqrt2 =2.*np.pi/(Dx*np.sqrt(2))

Vprime = np.zeros((gp_depth,gp_xrange))
Zprime = np.zeros((gp_depth,gp_xrange))
for iz in range(gp_depth):
    Vprime[iz,:] = np.sqrt(2)*U_scale*np.cos(k_over_sqrt2*xp_m)
    Zprime[iz,:] = -zeta_max*np.sin(k_over_sqrt2*xp_m)

#############Plot parameters###############
plot_slice=1
colormap='RdBu_r' 
aspect=0.6
ncases=1
nts=1

nxticks=2
txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)
nyticks=4
tylabels=np.arange(0,depth+1,depth/nyticks)
ticksy_loc=np.arange(0,gp_depth,(gp_depth-1.)/nyticks)

bu_min=0.
bu_max=2.
omega_min=-0.3
omega_max=0.3
wke_threshold = 1e-5
wke_levels=[1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
##########################################



#Create folder for plots if it doesn't exists                                                                                                                 
if not os.path.exists('plots/'+run+'/'+field_name[field-1]+'/'):
    os.makedirs('plots/'+run+'/'+field_name[field-1]+'/')

for ts in ts_plot:

    print "ts=",ts

    zeros_ts = (3-len(str(ts)))*'0'
    path_omega  = 'plots/'+run+'/omega/omega_t'+zeros_ts+str(ts)+'.dat'  
    path_wke    = 'plots/'+run+'/wn/wke_t'+zeros_ts+str(ts)+'.dat'                
    path_k      = 'plots/'+run+'/wn/k_t'+zeros_ts+str(ts)+'.dat'         
    path_l      = 'plots/'+run+'/wn/l_t'+zeros_ts+str(ts)+'.dat'         
    path_m      = 'plots/'+run+'/wn/m_t'+zeros_ts+str(ts)+'.dat'         

    #Load the data file                                                                                                                                                      
    if os.path.isfile(path_omega):

        omega = np.loadtxt(path_omega)
        k = np.loadtxt(path_k)
        l = np.loadtxt(path_l)
        m = np.loadtxt(path_m)
        wke = np.loadtxt(path_wke)

        Bu = (N2/f2)*(k*k+l*l)/(m*m)
        DS = Vprime*l/cor



        DS_norm  = DS/(zeta_max/2.)
        sig_norm = omega/(zeta_max/2.)
        Bu_norm  = 0.5*Bu/(zeta_max/2.)

        WKB_residual = np.absolute((omega-Zprime/2.-0.5*Bu - DS)/(zeta_max/2.))
#        WKB_residual = np.absolute(omega-Zprime/2.-0.5*Bu - DS)

        lambda_x = 2.*np.pi/k
        lambda_y = 2.*np.pi/l

        nan_indices = np.where(wke < 2.*wke_threshold)
        Bu[nan_indices] = np.nan
        DS[nan_indices] = np.nan
        Bu_norm[nan_indices] = np.nan
        DS_norm[nan_indices] = np.nan
        sig_norm[nan_indices] = np.nan
        omega[nan_indices] = np.nan
        lambda_x[nan_indices] = np.nan
        lambda_y[nan_indices] = np.nan
        WKB_residual[nan_indices] = np.nan


        slope = (omega-Zprime/2.)/(zeta_max/2.)

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
                
                if(field==1):
                    im = ax.imshow(sig_norm,cmap=colormap,aspect=aspect,vmin=-1.,vmax=1.)
                elif(field==2):
                    im = ax.imshow(Bu_norm,cmap=colormap,aspect=aspect,vmin=0.,vmax=2.)
                elif(field==3):                    
                    im = ax.imshow(DS_norm,cmap=colormap,aspect=aspect,vmin=-1.,vmax=1.)
                elif(field==4):                    
                    im = ax.imshow(WKB_residual,cmap=colormap,aspect=aspect,vmin=0.,vmax=1.)
                elif(field==5):
                    im = ax.imshow(slope,cmap=colormap,aspect=aspect,vmin=-0.25,vmax=0.25)                    
                elif(field==6):
                    #im = ax.imshow(Zprime/2.,cmap=colormap,aspect=aspect)                    
                    im = ax.imshow(omega,cmap=colormap,aspect=aspect)                    

#                im = ax.imshow(Bu,cmap=colormap,aspect=aspect,vmin=bu_min,vmax=bu_max)
#                im = ax.imshow(lambda_y,cmap=colormap,aspect=aspect,vmin=ly_min,vmax=ly_max)
#                im = ax.imshow(lambda_x,cmap=colormap,aspect=aspect,vmin=lx_min,vmax=lx_max)
#                im = ax.imshow(WKB_residual,cmap=colormap,aspect=aspect,vmin=0.,vmax=1.)
#                im = ax.imshow(Bu,cmap=colormap,aspect=aspect,vmin=0.,vmax=1.)
                cbar = ax.cax.colorbar(im)# ticks=[1e-8, 1e-2])                                                                                                          
                cbar = grid.cbar_axes[0].colorbar(im)
                
                time=ts*timestep
                time_title = '%.1f' % time
                ax.set_title(r''+field_title[field-1]+', $t =$ '+time_title+' inertial periods',fontsize=12)
                ax.text(-15, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
                ax.text(gp_del, gp_depth+20,r"$x'$ (km)",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)
                
                
                plt.show()                                                                                                                                                     
#                zeros_ts = (3-len(str(ts)))*'0'
#                plt.savefig('plots/'+run+'/'+field_name[field-1]+'/'+field_name[field-1]+'_t'+zeros_ts+str(ts)+'.png',bbox_inches='tight')

#make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/'+field_name[field-1]+'/*.png plots/'+run+'/'+field_name[field-1]+'/'+field_name[field-1]+'.gif'
#p = subprocess.Popen(make_gif, shell = True)
#os.waitpid(p.pid, 0)

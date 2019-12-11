import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os
import subprocess
import sys
from finds import find_resolution
from finds import find_scales
from finds import find_timestep

make_gif=1
show=0
plot_slice=1

leif_field=1
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'#double_gaussian/'#'leif/'



iy_transect = -int(256/4)
ix_offset   = int(iy_transect/2)
run = 'strain/N2_1e-5_y_n2-4_noadv_nodisp'#'real/ml_1024_uwz'
zoom=''#'_zoom'#''


location = scratch_location+folder+run
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)




depth = 1000  #m
#Range in x'
xpl = -35000#-350000#-35000#5000#
xpr = 35000#20000#
ts_list=np.arange(0,151,1)


field='dudz'#'dudz'#'dudz'

if(field=='dudz'):
    vmin = -0.0008
    vmax = 0.0008
elif(field=='wke'):
    vmin = 0.#0.001#-0.0008
    vmax = 0.001#0.0008
elif(field=='u'):
    vmin = -0.1#0.001#-0.0008
    vmax = 0.1#0.0008


timestep=0.1 #0.1 #Fraction of an inertial period between slices
hres=n1
vres=n3

colormap='RdBu_r' 
ncases=1
nts=1


#Translate location into grid points
dx=Dx/hres
dz=Dz/vres



x0 = int((Dx/2. + xpl*np.cos(np.deg2rad(45)))/dx) + ix_offset
x1 = int((Dx/2. + xpr*np.cos(np.deg2rad(45)))/dx) + ix_offset

print(x0,x1)

gp_del= x1-x0
gp_depth = int(vres*depth/Dz)  
aspect=0.5*gp_del/gp_depth

nxticks=2
txlabels=np.arange(xpl/1000,(xpr)/1000+1,(xpr-xpl)/(nxticks*1000))
ticksx_loc=np.arange(0,gp_del,(gp_del-1.)/nxticks)

nyticks=4
tylabels=np.arange(0,depth+1,depth/nyticks)
ticksy_loc=np.arange(0,gp_depth,(gp_depth-1.)/nyticks)




#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run+'/'+field+'/'):
    os.makedirs('plots/'+run+'/'+field+'/')

#u = np.zeros((gp_depth,x1-x0))
#v = np.zeros((gp_depth,x1-x0))
dudz = np.zeros((gp_depth,x1-x0))
dvdz = np.zeros((gp_depth,x1-x0))

for ts in ts_list:
    
    spaces_ts = (3-len(str(ts)))*' '
    path_lar  = scratch_location+folder+run+'/output/slicev1'+spaces_ts+str(ts)+'.dat'
    path_lai  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
    
    
    #Load the data file                                                                                                                  
    if os.path.isfile(path_lar):
    
        time=ts*timestep #time in inertial periods
        print('Time = ',time,' inertial periods.')

        if(leif_field==1):

            g_lar   = np.flipud(np.reshape(np.loadtxt(path_lar),(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                
            g_lai   = np.flipud(np.reshape(np.loadtxt(path_lai),(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                
        
        else:

#            g_lar   = np.rot90(np.reshape(np.loadtxt(path_lar),(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
#            g_lai   = np.rot90(np.reshape(np.loadtxt(path_lai),(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                

            g_lar   = np.flipud(np.fliplr(np.reshape(np.loadtxt(path_lar),(vres,hres))))        #Reshapes the array into a 2-d one  g(z,x,t)                                
            g_lai   = np.flipud(np.fliplr(np.reshape(np.loadtxt(path_lai),(vres,hres))))        #Reshapes the array into a 2-d one  g(z,x,t)                                


        u =  g_lar[0:gp_depth,x0:x1]*np.cos(time*2.*np.pi) + g_lai[0:gp_depth,x0:x1]*np.sin(time*2.*np.pi)
        v = -g_lar[0:gp_depth,x0:x1]*np.sin(time*2.*np.pi) + g_lai[0:gp_depth,x0:x1]*np.cos(time*2.*np.pi)
        wke = 0.5*(u*u+v*v)
        

        if(field=='dudz'):
            for iz in range(0,len(u[:,0])-1):
                dudz[iz,:] = (u[iz,:]-u[iz+1,:])/dz
                dvdz[iz,:] = (v[iz,:]-v[iz+1,:])/dz
    

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

            time_title = '%.1f' % time
            
            if(field=='dudz'):
                im = ax.imshow(dudz,cmap=colormap,aspect=aspect,vmin=vmin,vmax=vmax)
                ax.set_title(r'$du/dz$ (s$^{-1}$), $t =$ '+time_title+' inertial periods',fontsize=12)    
                cbar = grid.cbar_axes[0].colorbar(im,ticks=[vmin,vmin/2,0.,vmax/2.,vmax])
            elif(field=='wke'):
                im = ax.imshow(wke,cmap=colormap,aspect=aspect,vmin=vmin,vmax=vmax)
                ax.set_title(r'WKE (m/s)$^{2}$, $t =$ '+time_title+' inertial periods',fontsize=12)    
                cbar = grid.cbar_axes[0].colorbar(im,ticks=[0.,vmax/2.,vmax])                        
            elif(field=='u'):
                im = ax.imshow(u,cmap=colormap,aspect=aspect,vmin=vmin,vmax=vmax)
                ax.set_title(r'u (m/s), $t =$ '+time_title+' inertial periods',fontsize=12)    
                cbar = grid.cbar_axes[0].colorbar(im,ticks=[vmin,vmin/2,0.,vmax/2.,vmax])                        



#            cbar = ax.cax.colorbar(im)
#            cbar = grid.cbar_axes[0].colorbar(im)
        
            ax.text(-gp_del/5, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
            ax.text(gp_del/2, gp_depth+gp_depth/7,r"$x'$ (km)",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)

            if(show==1):
                plt.show()
            else:
                zeros_ts = (3-len(str(ts)))*'0'
                if(make_gif==1):
                    plt.savefig('plots/'+run+'/'+field+'/'+field+zoom+zeros_ts+str(ts)+'.png',bbox_inches='tight')
                else:
                    plt.savefig('plots/'+run+'/'+field+'/'+field+zoom+zeros_ts+str(ts)+'.eps',bbox_inches='tight')
    

if(make_gif==1):
    make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/'+field+'/'+field+zoom+'*.png plots/'+run+'/'+field+'/'+field+'.gif'
    p = subprocess.Popen(make_gif, shell = True)
    os.waitpid(p.pid, 0)

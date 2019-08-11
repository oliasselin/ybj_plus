import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

ts_list=[100]#np.arange(1,201,1)


field='dudz'#'dudz'#'dudz'

if(field=='dudz'):
    vmin = -0.001
    vmax = 0.001
else:
    vmin = 0.#0.001#-0.0008
    vmax = 0.001#0.0008

aspect=0.2

timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'#'leif/'
run = 'N2_1e5_a'#'attempt3_ro50'#'ml_100'

show=1
plot_slice=1
colormap='RdBu_r' 
ncases=1
nts=1

depth = 800  #m
Dz= 3000
gp_depth = int(vres*depth/Dz)
dz= Dz/vres

nyticks=4
tylabels=np.arange(0,depth+1,depth/nyticks)
ticksy_loc=np.arange(0,gp_depth,(gp_depth-1.)/nyticks)


Dx_km=100
Dx=Dx_km*1000
dx=Dx/hres

#Range in x'
xpl =10000#-35000#5000#
xpr =25000#20000#

#Corresponding grid point range
#x0 = int((Dx/2. - xpr*np.cos(np.deg2rad(45)))/dx)
#x1 = int((Dx/2. - xpl*np.cos(np.deg2rad(45)))/dx)
x0 = int((Dx/2. + xpl*np.cos(np.deg2rad(45)))/dx)
x1 = int((Dx/2. + xpr*np.cos(np.deg2rad(45)))/dx)
gp_del= x1-x0

print x0,x1

nxticks=3
txlabels=np.arange(xpl/1000,(xpr+1)/1000,(xpr-xpl)/(nxticks*1000))
ticksx_loc=np.arange(0,gp_del,(gp_del-1.)/nxticks)


#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run+'/'+field+'/'):
    os.makedirs('plots/'+run+'/'+field+'/')

u = np.zeros((gp_depth,x1-x0))
#v = np.zeros((gp_depth,x1-x0))
dudz = np.zeros((gp_depth,x1-x0))
#dvdz = np.zeros((gp_depth,x1-x0))

for ts in ts_list:
    
    spaces_ts = (3-len(str(ts)))*' '
    path_wke  = scratch_location+folder+run+'/output/slicev1'+spaces_ts+str(ts)+'.dat'
    path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
    path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'
    
    
    #Load the data file                                                                                                                  
    if os.path.isfile(path_wke):
    
        time=ts*timestep #time in inertial periods
        print 'Time = ',time,' inertial periods.'
        

        if(field == 'wke'):
            g_wke   = np.rot90(np.reshape(np.loadtxt(path_wke),(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
            wke = g_wke[0:gp_depth,x0:x1]          

        elif(field == 'dudz'):

            g_lar   = np.rot90(np.reshape(np.loadtxt(path_lar),(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
            g_lai   = np.rot90(np.reshape(np.loadtxt(path_lai),(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                

            u =  g_lar[0:gp_depth,x0:x1]*np.cos(time*2.*np.pi) + g_lai[0:gp_depth,x0:x1]*np.sin(time*2.*np.pi)
#            v = -g_lar[0:gp_depth,x0:x1]*np.sin(time*2.*np.pi) + g_lai[0:gp_depth,x0:x1]*np.cos(time*2.*np.pi)
                
            for iz in range(0,len(u[:,0])-1):
                dudz[iz,:] = (u[iz,:]-u[iz+1,:])/dz
#                dvdz[iz,:] = (v[iz,:]-v[iz+1,:])/dz
    

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


#            cbar = ax.cax.colorbar(im)
#            cbar = grid.cbar_axes[0].colorbar(im)
        
            ax.text(-3.5, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
            ax.text(gp_del/2, gp_depth+10,r"$x'$ (km)",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)

            if(show==1):
                plt.show()
            else:
                zeros_ts = (3-len(str(ts)))*'0'
                plt.savefig('plots/'+run+'/'+field+'/'+field+zeros_ts+str(ts)+'.eps',bbox_inches='tight')
    

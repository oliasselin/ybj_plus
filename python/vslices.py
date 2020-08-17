import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'N2_5e6'

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

xrange_km = 10
Dx_km=100
gp_del= int(hres*xrange_km/Dx_km) 
x0 = int(hres/2-gp_del) 
x1 = int(hres/2+gp_del) 

nxticks=4
txlabels=np.arange(-xrange_km,xrange_km+1,(2*xrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)








#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'):
    os.makedirs('plots/')

ts=100

g_wke = np.zeros((vres,hres))    
g_lar = np.zeros((vres,hres))    
g_lai = np.zeros((vres,hres))    

dudz = np.zeros((gp_depth,x1-x0))
dvdz = np.zeros((gp_depth,x1-x0))

spaces_ts = (3-len(str(ts)))*' '
path_wke  = scratch_location+folder+run+'/output/slicev1'+spaces_ts+str(ts)+'.dat'
path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'


    #Load the data file                                                                                                                  
if os.path.isfile(path_wke):
    
    time=ts*timestep #time in inertial periods
    print 'Time = ',time,' inertial periods.'

    f_wke = np.loadtxt(path_wke)                 #Loads the full file as a 1-dim array                                                      
    f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                      
    f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                      

    g_wke[:,:]   = np.rot90(np.reshape(f_wke,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
    g_lar[:,:]   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
    g_lai[:,:]   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                

    
    wke = g_wke[0:gp_depth,x0:x1]          
    lar = g_lar[0:gp_depth,x0:x1]          
    lai = g_lai[0:gp_depth,x0:x1]          
    
    u =  lar*np.cos(time*2.*np.pi) + lai*np.sin(time*2.*np.pi)
    v = -lar*np.sin(time*2.*np.pi) + lai*np.cos(time*2.*np.pi)


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
        
#        im = ax.imshow(dudz,cmap=colormap,aspect=0.4)
        im = ax.imshow(wke,cmap=colormap,aspect=1.)
#        im = ax.imshow(u,cmap=colormap)
#        im = ax.imshow(wke,cmap=colormap)
        
        cbar = ax.cax.colorbar(im)
        cbar = grid.cbar_axes[0].colorbar(im)
        



#        ax.set_title(r'$du/dz ($s$^{-1}$), $t =$ '+str(time)+' inertial periods',fontsize=12)    
        ax.set_title(r'WKE (m/s)$^{2}$, $t =$ '+str(time)+' inertial periods',fontsize=12)    
        ax.text(-7, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
        ax.text(gp_del, gp_depth+10,r'$x_{cs}$ (km)',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)


        plt.show()
#        plt.savefig('plots/'+run+'/dudz.eps',bbox_inches='tight')
#        plt.savefig('plots/'+run+'/wke.eps',bbox_inches='tight')
    

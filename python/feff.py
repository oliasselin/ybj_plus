import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'Dx_200'

plot_slice=1
colormap='RdBu_r' 
ncases=1
nts=1

aspect=0.4
depth = 3000#800  #m
Dz= 3000
gp_depth = int(vres*depth/Dz)
dz= Dz/vres

nyticks=4
tylabels=np.arange(0,depth+1,depth/nyticks)
ticksy_loc=np.arange(0,gp_depth,(gp_depth-1.)/nyticks)

xrange_km = 100#50#10
Dx_km=200#100
gp_del= int(hres*xrange_km/Dx_km) 
x0 = int(hres/2-gp_del) 
x1 = int(hres/2+gp_del) 
gp_xrange=x1-x0

nxticks=4
txlabels=np.arange(-xrange_km,xrange_km+1,(2*xrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)





#Time range (time should be in fractions of inertial period. Fraction = timestep var above)                                                                                                                   
ts_min=100
ts_max=151

#Declare vars 
u = np.zeros((gp_depth,gp_xrange,ts_max-ts_min))
v = np.zeros((gp_depth,gp_xrange,ts_max-ts_min))
time = np.zeros(ts_max-ts_min)
stairs = np.zeros((gp_depth,gp_xrange,ts_max-ts_min))



#Suboptimal: load vertical slice                                                                                                                                                                                     
for ts in range(ts_min,ts_max):

    its=ts-ts_min

    print its,ts

    spaces_ts = (3-len(str(ts)))*' '
    path_wke  = scratch_location+folder+run+'/output/slicev1'+spaces_ts+str(ts)+'.dat'
    path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
    path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'


    #Load the data file                                                                                                                                                                                    
    if os.path.isfile(path_wke):

        time[its]=ts*timestep #time in inertial periods         

        f_wke = np.loadtxt(path_wke)                 #Loads the full file as a 1-dim array                                                                                                                     
        f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                                                                                  
        f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                                                                                  

        g_wke   = np.rot90(np.reshape(f_wke,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                               
        g_lar   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                                      
        g_lai   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                                      

        lar = g_lar[0:gp_depth,x0:x1]                                                                                                                                                                        
        lai = g_lai[0:gp_depth,x0:x1]

        u[:,:,its] =  lar*np.cos(time[its]*2.*np.pi) + lai*np.sin(time[its]*2.*np.pi)
        v[:,:,its] = -lar*np.sin(time[its]*2.*np.pi) + lai*np.cos(time[its]*2.*np.pi)

        if(ts==(ts_min+int((ts_max-ts_min)/2))):
            wke = g_wke[0:gp_depth,x0:x1]                                                                                                                                                                    


        #Create a staircase function for both a and c points to make phi monotonuous                                                                                                                                   
        #Create a staircase function for both a and c points to make phi monotonuous                                                                                                                                    
        if its >= 1:
            for iz in range(len(u[:,0,0])):
                for ix in range(len(u[0,:,0])):

                    if np.arctan(v[iz,ix,its]/u[iz,ix,its]) > np.arctan(v[iz,ix,its-1]/u[iz,ix,its-1]):
                        stairs[iz,ix,its]=stairs[iz,ix,its-1]+np.pi
                    else:
                        stairs[iz,ix,its]=stairs[iz,ix,its-1]
        else:
            stairs[:,:,its]=0

#phi_final-phi_initial                                                                                                                                                                                                  
dphi = np.arctan(v[:,:,-1]/u[:,:,-1])-stairs[:,:,-1]  -  (np.arctan(v[:,:,0]/u[:,:,0])-stairs[:,:,0])
f_eff = - dphi/(2.*np.pi*(time[-1]-time[0]))


#Just calculate vorticity at the a and c points                                                                                                                                                        
path_vort  = scratch_location+folder+run+'/output/slice2v7 0.dat'
f_vort = np.loadtxt(path_vort)
g_vort = np.rot90(np.reshape(f_vort,(vres,hres)),k=2)

    

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
        im = ax.imshow(f_eff,cmap=colormap,aspect=aspect)
#        im = ax.imshow(u,cmap=colormap)
#        im = ax.imshow(wke,cmap=colormap)
        
        cbar = ax.cax.colorbar(im)
        cbar = grid.cbar_axes[0].colorbar(im)
        



#        ax.set_title(r'$du/dz ($s$^{-1}$), $t =$ '+str(time)+' inertial periods',fontsize=12)    
        ax.set_title(r'$f_{eff}/f$, $t =$ '+str(int(time[0]))+' - '+str(int(time[-1]))+'  inertial periods',fontsize=12)
#        ax.set_title(r'WKE (m/s)$^{2}$, $t =$ '+str(time)+' inertial periods',fontsize=12)    
        ax.text(-7, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
        ax.text(gp_del, gp_depth+10,r'$x_{cs}$ (km)',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)


#        plt.show()
#        plt.savefig('plots/'+run+'/dudz.eps',bbox_inches='tight')
#        plt.savefig('plots/'+run+'/wke.eps',bbox_inches='tight')
        plt.savefig('plots/'+run+'/feff_t'+str(int(time[0]))+'.eps',bbox_inches='tight')
    

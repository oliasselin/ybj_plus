import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

Dx = 100  #km (size of domain in the horizontal) 


scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'attempt3_ro50'

plot_slice=1
colormap='RdBu_r' 
ncases=1
nts=1

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run+'/'):
    os.makedirs('plots/'+run+'/')

leif_dist = 9 #km                                                                                                                                                                      
leif_dist_x = leif_dist*np.cos(np.deg2rad(45)) #Projection on the x-axis assuming 45 deg angle                                                                                   
gpdist = int(hres*leif_dist_x/Dx)

c_loc = hres/2-int(gpdist/2)          #location in x for the cyclone                                                                                                         
a_loc = hres/2+gpdist-int(gpdist/2)   #location in x for the anti-cyclone                                                                                                     




#Just calculate vorticity at the a and c points
path_vort  = scratch_location+folder+run+'/output/slice2htop7 0.dat'
f_vort = np.loadtxt(path_vort)
g_vort = np.rot90(np.reshape(f_vort,(hres,hres)),k=1) 

print np.max(g_vort)

path_u  = scratch_location+folder+run+'/output/slice2htop1 0.dat'
f_u = np.loadtxt(path_u)
g_u = np.rot90(np.reshape(f_u,(hres,hres)),k=0) 

path_v  = scratch_location+folder+run+'/output/slice2htop2 0.dat'
f_v = np.loadtxt(path_v)
g_v = np.rot90(np.reshape(f_v,(hres,hres)),k=0) 


#X = np.linspace(0,hres)
#Y = np.linspace(0,hres)

X, Y = np.meshgrid(np.arange(0,hres), np.arange(0,hres))

nticks=2
ticks_loc=np.arange(0,hres,(hres-1.)/nticks)
tlabels=np.arange(0,Dx+1,Dx/(nticks))

print ticks_loc
print tlabels


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
        
#        ax.get_xaxis().set_ticks([0 222])
#        ax.get_yaxis().set_ticks([0 222])
#        ax.get_xaxis().set_ticks([])
#        ax.get_yaxis().set_ticks([])
        
        ax.get_xaxis().set_ticks(ticks_loc)
        ax.set_xticklabels(tlabels)
        ax.get_xaxis().set_label('x (km)')

        ax.get_yaxis().set_ticks(ticks_loc)
        ax.set_yticklabels(tlabels[::-1])
        ax.get_yaxis().set_label('y (km)')

        ax.set_title(r'Surface $\zeta/f$')

        im = ax.imshow(g_vort,cmap=colormap)
        Q = ax.quiver(X[::16,::16],Y[::16,::16],g_u[::16,::16],g_v[::16,::16],color='k',pivot='mid', units='width')
        qk = ax.quiverkey(Q, 0.75, 0.92, 0.5, r'$50$ cm/s', labelpos='E',
                   coordinates='figure')

#        ax.scatter(c_loc, c_loc, color='m', s=20)
#        ax.scatter(a_loc, a_loc, color='y', s=20)


        cbar = ax.cax.colorbar(im)
        cbar = grid.cbar_axes[0].colorbar(im)
        
#        plt.show()
        plt.savefig('plots/'+run+'/vort.eps')


import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

Dx = 120  #km (size of domain in the horizontal) 

ts='100'

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'real_N2_1e5'#'real_N2_1e5'#'real_dg_ml'

plot_slice=1
field_to_plot='wke' #vort, u or v
show_plot=0
show_xyp=0
show_quivers=1
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
#g_vort = np.rot90(np.reshape(f_vort,(hres,hres)),k=0) 
g_vort = np.flipud(np.reshape(f_vort,(hres,hres))) 

path_u  = scratch_location+folder+run+'/output/slice2htop1 0.dat'
f_u = np.loadtxt(path_u)
#g_u = np.rot90(np.reshape(f_u,(hres,hres)),k=0) 
g_u = np.flipud(np.reshape(f_u,(hres,hres))) 

path_v  = scratch_location+folder+run+'/output/slice2htop2 0.dat'
f_v = np.loadtxt(path_v)
#g_v = np.rot90(np.reshape(f_v,(hres,hres)),k=0) 
g_v = np.flipud(np.reshape(f_v,(hres,hres))) 


path_p  = scratch_location+folder+run+'/output/slice2htop3 0.dat'
f_p = np.loadtxt(path_p)
g_p = np.flipud(np.reshape(f_p,(hres,hres))) 


if(field_to_plot=='wke'):
    path_lar = scratch_location+folder+run+'/output/slicehtop1'+ts+'.dat'
    f_lar = np.loadtxt(path_lar)
    g_lar = np.flipud(np.reshape(f_lar,(hres,hres)))

    path_lai = scratch_location+folder+run+'/output/slicehtop2'+ts+'.dat'
    f_lai = np.loadtxt(path_lai)
    g_lai = np.flipud(np.reshape(f_lai,(hres,hres)))

    wke = 0.5*(g_lar*g_lar + g_lai*g_lai)



X, Y = np.meshgrid(np.arange(0,hres), np.arange(0,hres))

nticks=2
ticks_loc=np.arange(0,hres,(hres-1.)/nticks)
tlabels=np.arange(-Dx/2,Dx/2+1,Dx/(nticks))


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



        if(field_to_plot=='vort'):
            ax.set_title(r'Surface $\zeta/f$')
            im = ax.imshow(g_vort,cmap=colormap)
        elif(field_to_plot=='u'):
            ax.set_title(r'Surface $u$ (m/s)')
            im = ax.imshow(g_u,cmap=colormap)
        elif(field_to_plot=='v'): 
            ax.set_title(r'Surface $v$ (m/s)')
            im = ax.imshow(g_v,cmap=colormap)
        elif(field_to_plot=='psi'): 
            ax.set_title(r'Surface $\psi$ (m$^2$/s)')
            im = ax.imshow(g_p,cmap=colormap)
        elif(field_to_plot=='wke'): 
            ax.set_title(r'Surface WKE (m/s)$^2$')
            im = ax.imshow(wke,cmap=colormap)

        if(show_quivers==1):
            Q = ax.quiver(X[::16,::16],Y[::16,::16],g_u[::16,::16],g_v[::16,::16],color='k',pivot='mid', units='width')
            qk = ax.quiverkey(Q, 0.75, 0.92, 0.5, r'$50$ cm/s', labelpos='E',
                              coordinates='figure')


        ax.text(-hres/10, hres/2,r'y (km)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
        ax.text(hres/2, hres+hres/10,r"$x$ (km)",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)

        if(show_xyp==1):
            color='green'

            ax.arrow(hres/2,hres/2,hres/10,hres/10,width=0.5,head_width=5.,color=color)
            ax.arrow(hres/2,hres/2,hres/10,-hres/10,width=0.5,head_width=5.,color=color)

            ax.text(hres/2+hres/25,hres/2+hres/8,r"$x'$",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=18,color=color)
            ax.text(hres/2+hres/25,hres/2-hres/10,r"$y'$",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=18,color=color)            

        cbar = ax.cax.colorbar(im)
        cbar = grid.cbar_axes[0].colorbar(im)

        if(show_plot==1):
            plt.show()
        else:
            plt.savefig('plots/'+run+'/'+field_to_plot+'.eps')



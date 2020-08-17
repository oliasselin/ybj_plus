import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogLocator
import numpy as np
#from numpy import ma
import os
import subprocess
import sys
from finds import find_resolution
from finds import find_scales

field_to_plot='refraction' #vort, u or v
vloc='htop'
term_max = 5e-7
make_gif=1
plot_slice=1
show_plot=0
show_xyp=0
show_quivers=0
show_quivers_key=0
show_contours=0
log_plot=0

leif_field=1
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'#double_gaussian/'#'leif/'                                                                                                                                              
run = 'strain/N2_1e-5_y_n2-4_nodisp'#'attempt3_ro50'#'ml_100'                                                                                                                                   

location = scratch_location+folder+run
n1,n2,n3 = find_resolution(location)
Dx,Dz,L_scale,H_scale,U_scale,h_thermo = find_scales(location)

timestep=0.1   #Number of inertial periods between frames
ts_list=np.arange(0,131,10)                                                                                                                                                     
wke_min=0.
wke_max=0.005

vort_levels=[-0.2, -0.1, 0. ,0.1, 0.2]
wke_levels=np.logspace(start=-1,stop=1,num=3,base=10)#[1e-3,1e-2, 1e-1, 1e0, 1e1]

#Automatic parameter setting
colormap='RdBu_r' 
ncases=1
nts=1
hres=n1
vres=n3

nticks=2
ticks_loc=np.arange(0,hres,(hres-1.)/nticks)
tlabels=np.arange(-int(Dx/2/1000),int(Dx/2/1000+1),int(Dx/(nticks)/1000))


#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run+'/'+field_to_plot+'/'):
    os.makedirs('plots/'+run+'/'+field_to_plot+'/')

if(show_quivers==1):

    path_u  = scratch_location+folder+run+'/output/slice2htop1 0.dat'
    g_u = np.flipud(np.reshape(np.loadtxt(path_u),(hres,hres))) 

    path_v  = scratch_location+folder+run+'/output/slice2htop2 0.dat'
    g_v = np.flipud(np.reshape(np.loadtxt(path_v),(hres,hres))) 

    if(leif_field!=1):
        g_u = np.fliplr(g_u)
        g_v = np.fliplr(g_v)

    X, Y = np.meshgrid(np.arange(0,hres), np.arange(0,hres))


if(show_contours==1):

    path_vort  = scratch_location+folder+run+'/output/slice2htop7 0.dat'
    g_vort = np.flipud(np.reshape(np.loadtxt(path_vort),(hres,hres)))

    if(leif_field!=1):
        g_vort = np.fliplr(g_vort)

for ts in ts_list:

    time=ts*timestep #time in inertial periods                                                                                                                                      
    print('Time = ',time,' inertial periods.')


    if(field_to_plot=='wke'):
        spaces_ts = (3-len(str(ts)))*' '
        path_lar  = scratch_location+folder+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
        path_lai  = scratch_location+folder+run+'/output/slicehtop2'+spaces_ts+str(ts)+'.dat'

        g_lar = np.flipud(np.reshape(np.loadtxt(path_lar),(hres,hres)))
        g_lai = np.flipud(np.reshape(np.loadtxt(path_lai),(hres,hres)))

        if(leif_field!=1):
            g_lar = np.fliplr(g_lar)
            g_lai = np.fliplr(g_lai)

        field = 0.5*(g_lar*g_lar + g_lai*g_lai)


    if(field_to_plot=='wke_vave'):
        spaces_ts = (3-len(str(ts)))*' '
        path_wvave  = scratch_location+folder+run+'/output/WE_vave'+spaces_ts+str(ts)+'.dat'

        field = np.flipud(np.reshape(np.loadtxt(path_wvave),(hres,hres)))

        ave_energy = np.average(field)
        field = field/ave_energy

        if(leif_field!=1):
            field = np.fliplr(field)


    if(field_to_plot=='tendency' or field_to_plot=='advection' or field_to_plot=='refraction' or field_to_plot=='dispersion'):

        if(field_to_plot=='tendency'):
            field_no = 1
        if(field_to_plot=='advection'):
            field_no = 3
        if(field_to_plot=='refraction'):
            field_no = 5
        if(field_to_plot=='dispersion'):
            field_no = 7

        spaces_ts = (3-len(str(ts)))*' '
        path_re  = scratch_location+folder+run+'/output/slice'+vloc+'3'+str(field_no)+spaces_ts+str(ts)+'.dat'
        path_im  = scratch_location+folder+run+'/output/slice'+vloc+'3'+str(field_no+1)+spaces_ts+str(ts)+'.dat'

        field_re = np.flipud(np.reshape(np.loadtxt(path_re),(hres,hres)))
        field_im = np.flipud(np.reshape(np.loadtxt(path_im),(hres,hres)))

        field = np.sqrt( np.square(field_re)+np.square(field_im) )

        if(leif_field!=1):
            field = np.fliplr(field)




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
        
            ax.get_xaxis().set_ticks(ticks_loc)
            ax.set_xticklabels(tlabels)
            ax.get_xaxis().set_label('x (km)')
            
            ax.get_yaxis().set_ticks(ticks_loc)
            ax.set_yticklabels(tlabels[::-1])
            ax.get_yaxis().set_label('y (km)')

            time_title = '%.1f' % time

            if(field_to_plot=='wke'): 
                ax.set_title(r'Surface WKE (m/s)$^2$, $t =$ '+time_title+' IP',fontsize=12)
                im = ax.imshow(field,cmap=colormap,vmin=wke_min,vmax=wke_max)

            if(field_to_plot=='wke_vave'):
                ax.set_title(r'Vert.-ave WE anomaly, $t =$ '+time_title+' IP',fontsize=12)
                if(log_plot==1):
                    im = ax.imshow(field,cmap=colormap,norm=LogNorm(vmin=wke_levels[0], vmax=wke_levels[-1]))#vmin=0.)#,vmax=10.)
                else:
                    im = ax.imshow(field,cmap=colormap,vmin=0.,vmax=10.) 
                
            if(field_to_plot=='tendency'):
                ax.set_title(r'Tendency (m/s$^2$), $t =$ '+time_title+' IP',fontsize=12)
                im = ax.imshow(field,cmap=colormap,vmin=0,vmax=term_max)
            if(field_to_plot=='advection'):
                ax.set_title(r'Advection (m/s$^2$), $t =$ '+time_title+' IP',fontsize=12)
                im = ax.imshow(field,cmap=colormap,vmin=0,vmax=term_max)#,vmin=wke_min,vmax=wke_max)
            if(field_to_plot=='refraction'):
                ax.set_title(r'Refraction (m/s$^2$), $t =$ '+time_title+' IP',fontsize=12)
                im = ax.imshow(field,cmap=colormap,vmin=0,vmax=term_max)#,vmin=wke_min,vmax=wke_max)
            if(field_to_plot=='dispersion'):
                ax.set_title(r'Dispersion (m/s$^2$), $t =$ '+time_title+' IP',fontsize=12)
                im = ax.imshow(field,cmap=colormap,vmin=0,vmax=term_max)#,vmin=wke_min,vmax=wke_max)

            if(show_quivers==1):
                Q = ax.quiver(X[::16,::16],Y[::16,::16],g_u[::16,::16],g_v[::16,::16],color='k',pivot='mid', units='width')
                if(show_quivers_key==1):
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
  
            if(log_plot==1):
                cbar = ax.cax.colorbar(im,ticks=wke_levels)
            else:
                cbar = ax.cax.colorbar(im)
#            cbar = ax.cax.colorbar(im)#,ticks=wke_levels)
#            cbar.ax.yaxis.set_major_locator(LogLocator())  # <- Why? See above.
#            cbar.set_ticks(cbar.ax.yaxis.get_major_locator().tick_values(z.min(), z.max())
            
#            cbar.ax.set_yticklabels(['yo','sup','vachier','salut'])
#            cbar = ax.cax.colorbar(im)

#            cbar = grid.cbar_axes[0].colorbar(im)

            if(show_contours==1):
                im=ax.contour(g_vort,levels=vort_levels,colors='k')#,colors='k')   


            if(show_plot==1):
                plt.show()
            else:
                zeros_ts = (3-len(str(ts)))*'0'
                if(make_gif==1):
                    plt.savefig('plots/'+run+'/'+field_to_plot+'/'+field_to_plot+zeros_ts+str(ts)+'.png')
                    print('plots/'+run+'/'+field_to_plot+zeros_ts+str(ts)+'.png')
                else:
                    plt.savefig('plots/'+run+'/'+field_to_plot+'/'+field_to_plot+zeros_ts+str(ts)+'.eps')


if(make_gif==1):
    make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/'+field_to_plot+'/*.png plots/'+run+'/'+field_to_plot+'/'+field_to_plot+'.gif'
    p = subprocess.Popen(make_gif, shell = True)
    os.waitpid(p.pid, 0)


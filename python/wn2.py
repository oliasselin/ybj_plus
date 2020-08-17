import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os
import subprocess
import sys

timestep=0.1 #0.1 #Fraction of an inertial period between slices
cor=1.24e-4


hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'

iy_transect = -int(hres/4)
ix_offset = int(iy_transect/2)
run =  'strain/N2_1e-5_y_n2-4_noadv_nodisp'#'N2_2e5_a'


leif_field=1


plot_slice=0
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
x0 = int(hres/2-gp_del) + ix_offset
x1 = int(hres/2+gp_del) + ix_offset
gp_xrange=x1-x0
#dx=Dx_km*np.sqrt(2)/hres#Dx_km*np.cos(np.deg2rad(45))/hres

nxticks=2
#txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)
###

lv_min=-300
lv_max=300
lh_min=-30
lh_max=30
wke_threshold = 1e-5

#For plots only
wke_min=1e-5
wke_max=0.01

#Time range (time should be in fractions of inertial period. Fraction = timestep var above)                                                                                                                   
ts_min=0
ts_max=161


#Create folder for plots if it doesn't exists                                                                                                                                            
if not os.path.exists('plots/'+run+'/wn/'):
    os.makedirs('plots/'+run+'/wn/')
#Suboptimal: load vertical slice                                                                                                                              

  
dlardz = np.zeros((gp_depth,x1-x0))
dlaidz = np.zeros((gp_depth,x1-x0))
dlardx = np.zeros((gp_depth,x1-x0))
dlaidx = np.zeros((gp_depth,x1-x0))
  
for ts in range(ts_min,ts_max):

    
    print("ts=",ts)

    spaces_ts = (3-len(str(ts)))*' '
    path_lar  = scratch_location+folder+run+'/output/slicev1'+spaces_ts+str(ts)+'.dat'
    path_lai  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
    path_lart  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'
    path_lait  = scratch_location+folder+run+'/output/slicev4'+spaces_ts+str(ts)+'.dat'


    path_larx  = scratch_location+folder+run+'/output/slicev5'+spaces_ts+str(ts)+'.dat'
    path_laix  = scratch_location+folder+run+'/output/slicev6'+spaces_ts+str(ts)+'.dat'
    path_lary  = scratch_location+folder+run+'/output/slicev7'+spaces_ts+str(ts)+'.dat'
    path_laiy  = scratch_location+folder+run+'/output/slicev8'+spaces_ts+str(ts)+'.dat'


    #Load the data file                                                                                                                                                      
    if os.path.isfile(path_lar):

        time=ts*timestep #time in inertial periods         

        f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                                             
        f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                                                   

        f_lart = np.loadtxt(path_lart)                 #Loads the full file as a 1-dim array                                                                             
        f_lait = np.loadtxt(path_lait)                 #Loads the full file as a 1-dim array
        f_larx = np.loadtxt(path_larx)                 #Loads the full file as a 1-dim array                                                                             
        f_laix = np.loadtxt(path_laix)                 #Loads the full file as a 1-dim array
        f_lary = np.loadtxt(path_lary)                 #Loads the full file as a 1-dim array                                                                             
        f_laiy = np.loadtxt(path_laiy)                 #Loads the full file as a 1-dim array                                                                                   

        if(leif_field==1):

            g_lar   = np.flipud(np.reshape(f_lar,(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                                       
            g_lai   = np.flipud(np.reshape(f_lai,(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                                        

            g_lart  = np.flipud(np.reshape(f_lart,(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                                                  
            g_lait  = np.flipud(np.reshape(f_lait,(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                                                  
            g_larx  = np.flipud(np.reshape(f_larx,(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                                                 
            g_laix  = np.flipud(np.reshape(f_laix,(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                                                  
            g_lary  = np.flipud(np.reshape(f_lary,(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)                                                                 
            g_laiy  = np.flipud(np.reshape(f_laiy,(vres,hres)))        #Reshapes the array into a 2-d one  g(z,x,t)

        else:


            g_lar   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                   
            g_lai   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                         
        
            g_lart   = np.rot90(np.reshape(f_lart,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                
            g_lait   = np.rot90(np.reshape(f_lait,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                         
            g_larx   = np.rot90(np.reshape(f_larx,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                 
            g_laix   = np.rot90(np.reshape(f_laix,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                         
            g_lary   = np.rot90(np.reshape(f_lary,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                                 
            g_laiy   = np.rot90(np.reshape(f_laiy,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                                         
            
        lar = g_lar[0:gp_depth,x0:x1]                                                                                                                                
        lai = g_lai[0:gp_depth,x0:x1]

        lart = g_lart[0:gp_depth,x0:x1]                                                                                                                                
        lait = g_lait[0:gp_depth,x0:x1]
        larx = g_larx[0:gp_depth,x0:x1]                                                                                                                                
        laix = g_laix[0:gp_depth,x0:x1]
        lary = g_lary[0:gp_depth,x0:x1]                                                                                                                                
        laiy = g_laiy[0:gp_depth,x0:x1]

        for iz in range(1,len(lar[:,0])-1):
            dlardz[iz,:] = (lar[iz-1,:]-lar[iz+1,:])/(2.*dz)
            dlaidz[iz,:] = (lai[iz-1,:]-lai[iz+1,:])/(2.*dz)

#        for ix in range(1,len(lar[0,:])-1):
#            dlardx[:,ix] = (lar[:,ix+1]-lar[:,ix-1])/(2.*dx)
#            dlaidx[:,ix] = (lai[:,ix+1]-lai[:,ix-1])/(2.*dx)

        #Calculate the Eulerian frequency
        sig = (lart*lai-lait*lar)/(lar*lar+lai*lai)
        sig = sig/cor

        #Calculate the local wavenumbers in m^{-1}. Primed coordinate system: k' and l' are the cross- and along-stream wavenumbers.
        m = (dlaidz*lar - dlardz*lai)/(lar*lar+lai*lai)
        lambda_z = 2.*np.pi/m

        lambda_z[0,:] = 0.
        lambda_z[-1,:] = 0.

        #Using the spectral accurate LAx and LAy
        if(leif_field==1):
            k = ((laix-laiy)*lar + (lary-larx)*lai )/(np.sqrt(2)*(lar*lar+lai*lai))
        else:
            k = -((laix+laiy)*lar - (larx+lary)*lai )/(np.sqrt(2)*(lar*lar+lai*lai))
        lambda_x = 2.*np.pi/(k*1000)

        if(leif_field==1):
            l = ((laix+laiy)*lar - (larx+lary)*lai )/(np.sqrt(2)*(lar*lar+lai*lai))
        else:
            l = ((laix-laiy)*lar - (larx-lary)*lai )/(np.sqrt(2)*(lar*lar+lai*lai))
        lambda_y = 2.*np.pi/(l*1000)

        wke = 0.5*(lar*lar+lai*lai)


        zeros_ts = (3-len(str(ts)))*'0'
        np.savetxt('plots/'+run+'/wn/sig_t'+zeros_ts+str(ts)+'.dat',sig)
        np.savetxt('plots/'+run+'/wn/k_t'+zeros_ts+str(ts)+'.dat',k)
        np.savetxt('plots/'+run+'/wn/l_t'+zeros_ts+str(ts)+'.dat',l)
        np.savetxt('plots/'+run+'/wn/m_t'+zeros_ts+str(ts)+'.dat',m)
        np.savetxt('plots/'+run+'/wn/wke_t'+zeros_ts+str(ts)+'.dat',wke)

 


        if(plot_slice==1):

            nan_indices = np.where(lar*lar + lai*lai < 2.*wke_threshold)
            sig[nan_indices] = 0.
            lambda_z[nan_indices] = 0.
            lambda_x[nan_indices] = 0.
            lambda_y[nan_indices] = 0.

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

#               im = ax.imshow(wke,cmap=colormap,norm=colors.LogNorm(vmin=wke_min, vmax=wke_max),aspect=aspect)
                im = ax.imshow(lambda_x,cmap=colormap,vmin=lh_min,vmax=lh_max,aspect=aspect)
#                im = ax.imshow(lambda_z,cmap=colormap,vmin=lv_min,vmax=lv_max,aspect=aspect)
#                im = ax.imshow(sig,cmap=colormap,aspect=aspect)#,vmin=sig_min,vmax=sig_max)
                
#                im = ax.imshow(k,cmap=colormap,vmin=-0.0001,vmax=0.0001,aspect=aspect)
                #im = ax.imshow(wke,cmap=colormap,aspect=aspect)


                cbar = ax.cax.colorbar(im)# ticks=[1e-8, 1e-2])
                cbar = grid.cbar_axes[0].colorbar(im)

                time_title = '%.1f' % time        
#                ax.set_title(r'$\sigma/f$ (km), $t =$ '+time_title+' inertial periods',fontsize=12)
                ax.set_title(r'$\lambda_x$ (km), $t =$ '+time_title+' inertial periods',fontsize=12)
                ax.text(-15, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
                ax.text(gp_del, gp_depth+20,r'$x_{cs}$ (km)',rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)


                plt.show()
#                zeros_ts = (3-len(str(ts)))*'0'
#                plt.savefig('plots/'+run+'/lh/lh_t'+zeros_ts+str(ts)+'.png',bbox_inches='tight')

#if(plot_slice==1):
#    make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/lh/*.png plots/'+run+'/lh/lh.gif'
#    p = subprocess.Popen(make_gif, shell = True)
#    os.waitpid(p.pid, 0)

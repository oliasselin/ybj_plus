import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os
import subprocess
import sys

ts_plot=[100]#np.arange(5,101,5)#[30]#np.arange(0,180,1)#[66]#np.arange(0,181,10)#np.arange(141,201,1)#[50,100,150,200]

field = 5
trace = 0
wke_contours=1
dzdx  = 'km'
plot_slice=1
make_gif=0
show=1

field_title = ['$\sigma/(\zeta_{max}/2)$','$Bu$','Doppler shift $/(\zeta_{max}/2)$','$R_{WKB}$','$|\sigma-\zeta/2|/(\zeta_{max}/2)$','$\sigma/f$','$\zeta/\zeta_{max}$','-k/m','$Bu_p$','|$Bu-Bu_p$|','$|(Bu-Bu_p)(1-$min$($WKB$_{res},1))|$']
field_name  = ['sig','bu','ds','res','refractive_phase','omega','zeta','km','bu_pred','bu_error','bu_error_weighted']

if(trace==1):
    trace_name=dzdx
else:
    trace_name=''


#Run specification##############################################
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'N2_1e5_a' #'N2_1e5_dla'

N2=1e-5
cor = 1.24e-4
f2=np.power(cor,2)
zeta_max=2*np.pi/(100000*1.24e-4)  #Simply = k/f
grad_zeta=2.7e-9
U_scale = 0.5

hres=256
vres=256
timestep=0.1 #0.1 #Fraction of an inertial period between slices
################################################################



depth = 100                   #Depth of interest (m)
Dz= 3000                       #Full depth of the simulation (m)
gp_depth = int(vres*depth/Dz)  #Number of gp in the vertical
dz= Dz/vres                    #dz (m)


xyrange_km = 35#np.sqrt(2)*25                    #Half-range of x' (km)
xrange_km  = xyrange_km*np.cos(np.deg2rad(45))   #Half-range of x  (km) -- projection assuming 45 degrees angle
Dx_km=100                                        #Domain size (km)
Dx=Dx_km*1000                                    #Domain size (m)
gp_del= int(round(hres*xrange_km/Dx_km))         #Half the number of gp in x
x0 = int(hres/2-gp_del)                          #ID of the first gp   
x1 = int(hres/2+gp_del)                          #ID of the first gp   
gp_xrange=x1-x0                                  #Number of gps in x
dx = xyrange_km*1000.*2./(1.*gp_xrange)          #dx' (m)  

gpx = np.arange(0,gp_xrange,1)                   #array of gp in x'   
gpz = np.arange(0,gp_depth,1)                    #array of gp in z   

xp_m = -xyrange_km*1000 + dx*gpx                 #array of gp in m

k_over_sqrt2 =2.*np.pi/(Dx*np.sqrt(2))           #k/sqrt(2)

Vprime = np.zeros((gp_depth,gp_xrange))          #V'(z,x')    --- along-stream velocity)
Zprime = np.zeros((gp_depth,gp_xrange))          #Zeta'(z,x') --- vertical vorticity 
GradZp = np.zeros((gp_depth,gp_xrange))          #d Zeta'(z,x')/dx' --- gradient of vertical vorticity at y'=0                                                                          
bigD   = np.zeros((gp_depth,gp_xrange))          #Depth vector (>0)
for iz in range(gp_depth):
    Vprime[iz,:] = np.sqrt(2)*U_scale*np.cos(k_over_sqrt2*xp_m)
    Zprime[iz,:] = -zeta_max*np.sin(k_over_sqrt2*xp_m)
    GradZp[iz,:] = -grad_zeta*np.cos(k_over_sqrt2*xp_m)
    bigD[iz,:]   = dz*gpz[iz]
    
alpha = -0.5*GradZp   #k=alpha*t  (>0 in the neighborhood of x'=0)

#############Plot parameters###############
colormap='RdBu_r' 
aspect=1.*gp_del/gp_depth #0.6
ncases=1
nts=1

nxticks=2
txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)
nyticks=4
tylabels=np.arange(0,depth+1,depth/nyticks)
ticksy_loc=np.arange(0,gp_depth,(gp_depth-1.)/nyticks)

wke_threshold = 1e-10#1e-5#1e-4
wke_levels=[1e-5, 1e-4, 1e-3, 1e-2] #[1e-4, 1e-3, 1e-2]#[1e-5, 1e-4, 1e-3, 1e-2]
##########################################



####### Tracing starting points ##########
xstart_list=[int(1*gp_del/4),int(gp_del/2),int(3*gp_del/4),gp_del,int(5*gp_del/4),int(3*gp_del/2),int(7*gp_del/4)]#[(int(gp_xrange/2)-20)]
zstart=1


#Create folder for plots if it doesn't exists                                                                                                                 
if not os.path.exists('plots/'+run+'/'+field_name[field-1]+trace_name+'/'):
    os.makedirs('plots/'+run+'/'+field_name[field-1]+trace_name+'/')

for ts in ts_plot:

    print "ts=",ts

    zeros_ts = (3-len(str(ts)))*'0'
    path_sig    = 'plots/'+run+'/wn/sig_t'+zeros_ts+str(ts)+'.dat'  
    path_wke    = 'plots/'+run+'/wn/wke_t'+zeros_ts+str(ts)+'.dat'                
    path_k      = 'plots/'+run+'/wn/k_t'+zeros_ts+str(ts)+'.dat'         
    path_l      = 'plots/'+run+'/wn/l_t'+zeros_ts+str(ts)+'.dat'         
    path_m      = 'plots/'+run+'/wn/m_t'+zeros_ts+str(ts)+'.dat'         

    #Load the data file                                                                                                                                                      
    if os.path.isfile(path_sig):

        #Load frequency and wavenumbers from WKB and also WKE
        sig = np.loadtxt(path_sig)
        k = np.loadtxt(path_k)
        l = np.loadtxt(path_l)
        m = np.loadtxt(path_m)
        wke = np.loadtxt(path_wke)

        sig = sig[:gp_depth,:gp_xrange]
        k   =   k[:gp_depth,:gp_xrange]
        l   =   l[:gp_depth,:gp_xrange]
        m   =   m[:gp_depth,:gp_xrange]
        wke = wke[:gp_depth,:gp_xrange]

        #Let's do some ray tracing
        if(trace==1):


            #Trajectories in grid points
            ix_traj=np.zeros((len(xstart_list),gp_xrange))            
            iz_traj=np.zeros((len(xstart_list),gp_xrange))            

            #Compute trajectories                                                                                                                                                            
            for traj,xstart in enumerate(xstart_list):
            
                iz_traj[traj,0]=zstart
                z = - dz*(iz_traj[traj,0]+0.5)

                for it in range(0,gp_xrange-xstart):
                
                    ix_traj[traj,it]=xstart+it

                #Calculate slope
#                arg=(f2/N2)*(2.*sig[iz_traj[it],ix_traj[it]]-Zprime[iz_traj[it],ix_traj[it]])     #YBJ formulation

                    if(dzdx=='leif'):
                        sig2  = np.power(1+sig[iz_traj[traj,it],ix_traj[traj,it]],2)
                        feff2 = np.power(1+0.5*Zprime[iz_traj[traj,it],ix_traj[traj,it]],2)
                        arg=(f2/N2)*(sig2-feff2)
                    elif(dzdx=='ybj'):
                        arg=(f2/N2)*(2.*sig[iz_traj[traj,it],ix_traj[traj,it]]-Zprime[iz_traj[traj,it],ix_traj[traj,it]])
                    elif(dzdx=='km'): #k/m# 
                        arg=k[iz_traj[traj,it],ix_traj[traj,it]]/m[iz_traj[traj,it],ix_traj[traj,it]]

    
                    if(arg>=0):
                        #All is good! Proceed to calculate the slope
                        if(dzdx=='km'):
                            slope_ybj=-arg
                        else:
                            slope_ybj=-np.sqrt(arg)
                    elif(arg<0):
                        print "sig < f_eff. Stoping trajectory"
                        iz_traj[traj,it:]=iz_traj[traj,it]
                        ix_traj[traj,it:]=ix_traj[traj,it]
                        break
                    elif(np.isnan(arg)):
                        print "sig is NaN. Stoping trajectory"
                        iz_traj[traj,it:]=iz_traj[traj,it]
                        ix_traj[traj,it:]=ix_traj[traj,it]
                        break

                    z = z + slope_ybj*dx

                    #Update iz and stop if bottom is reached
                    iz_traj[traj,it+1]=-int(0.5+z/dz)
                    if(iz_traj[traj,it+1]>gp_depth-1):
                        print "Reached the bottom"
                        iz_traj[traj,it:]=iz_traj[traj,it]
                        ix_traj[traj,it:]=ix_traj[traj,it]
                        break

                    #Stop if rightmost limit is reached
                    if(it==(gp_xrange-xstart-1)):
                        print "Reached the anticyclonic core"
                        iz_traj[traj,it:]=iz_traj[traj,it]
                        ix_traj[traj,it:]=ix_traj[traj,it]
                        break

        #Remove regions where WKE is too weak.
        nan_indices = np.where(wke < wke_threshold)
        sig[nan_indices] = np.nan
        k[nan_indices] = np.nan
        l[nan_indices] = np.nan
        m[nan_indices] = np.nan

        #Calculate the WKB-based Burger number and Doppler shift (y'=0)
        Bu = (N2/f2)*(k*k+l*l)/(m*m)
        DS = Vprime*l/cor


        #Calculate the prediction to Bu
        Bu_pred = (N2/f2)*alpha*alpha*np.power(N2*alpha*alpha/(3.*cor*bigD),-2./3.)
        Bu_error= np.abs(Bu-Bu_pred)


        #Normalize by zeta_max/2
        DS_norm  = DS/(zeta_max/2.)
        sig_norm = sig/(zeta_max/2.)
        Bu_norm  = 0.5*Bu/(zeta_max/2.)

        #Other quantities of interest
        WKB_residual = np.absolute((sig-Zprime/2.-0.5*Bu - DS)/(zeta_max/2.))
        WKBness      = (1.-np.minimum(WKB_residual,WKB_residual*0.+1.))

        print WKB_residual.shape,WKBness.shape

        slope = (sig-Zprime/2.)/(zeta_max/2.)
#        slope = np.real(np.sqrt(   (f2/N2)*(sig*sig-np.power(Zprime/2.,2)) ))

        zeta = Zprime/zeta_max
        zeta[nan_indices] = np.nan


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
                    im = ax.imshow(Bu,cmap=colormap,aspect=aspect,vmin=0.,vmax=1.)
                elif(field==3):                    
                    im = ax.imshow(DS_norm,cmap=colormap,aspect=aspect,vmin=-1.,vmax=1.)
                elif(field==4):                    
                    im = ax.imshow(WKB_residual,cmap=colormap,aspect=aspect,vmin=0.,vmax=1.)
                elif(field==5):
#                    im = ax.imshow(slope,cmap=colormap,aspect=aspect,vmin=-1.,vmax=1.)                    
                    im = ax.imshow(np.abs(slope),cmap=colormap,aspect=aspect,vmin=0.,vmax=1.0)                    
                elif(field==6):
                    im = ax.imshow(sig,cmap=colormap,aspect=aspect,vmin=-0.3,vmax=0.3)
                    #im = ax.imshow(Zprime/2.,cmap=colormap,aspect=aspect)                    
                elif(field==7):
                    im = ax.imshow(zeta,cmap=colormap,aspect=aspect,vmin=-1.,vmax=1.)
                elif(field==8):
                    im = ax.imshow(-k/m,cmap=colormap,aspect=aspect,vmin=-0.04,vmax=0.0)
                elif(field==9):
                    im = ax.imshow(Bu_pred,cmap=colormap,aspect=aspect,vmin=0.,vmax=1.0)
                elif(field==10):
                    im = ax.imshow(Bu_error,cmap=colormap,aspect=aspect,vmin=0.,vmax=1.0)
                elif(field==11):
                    im = ax.imshow(np.abs(Bu_error*WKBness),cmap=colormap,aspect=aspect,vmin=0.,vmax=1.0)




                if(wke_contours==1):
                    ax.contour(wke,levels=wke_levels,colors='k')#, levels, colors='k', origin='upper', extent=extent)                                              


                cbar = ax.cax.colorbar(im)# ticks=[1e-8, 1e-2])                                                                                                          
                cbar = grid.cbar_axes[0].colorbar(im)
                
                time=ts*timestep
                time_title = '%.1f' % time
                ax.set_title(r''+field_title[field-1]+', $t =$ '+time_title+' inertial periods',fontsize=12)
                ax.text(-15, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
                ax.text(gp_del, gp_depth+gp_depth/7,r"$x'$ (km)",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)
                
                
#                for traj,xstart in enumerate(xstart_list):

                if(trace==1):
                    for traj,xstart in enumerate(xstart_list):
                        ax.plot(ix_traj[traj,:],iz_traj[traj,:],'-k',linewidth=2.5)

                if(show==1):
                    plt.show()
                else:
                    zeros_ts = (3-len(str(ts)))*'0'
                    if(make_gif==1):
                        plt.savefig('plots/'+run+'/'+field_name[field-1]+trace_name+'/'+field_name[field-1]+trace_name+'_t'+zeros_ts+str(ts)+'.png',bbox_inches='tight')
                    else:
                        plt.savefig('plots/'+run+'/'+field_name[field-1]+trace_name+'/'+field_name[field-1]+trace_name+'_t'+zeros_ts+str(ts)+'.eps',bbox_inches='tight')
            
                    
if(make_gif==1):
    make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/'+field_name[field-1]+trace_name+'/*.png plots/'+run+'/'+field_name[field-1]+trace_name+'/'+field_name[field-1]+trace_name+'.gif'
    p = subprocess.Popen(make_gif, shell = True)
    os.waitpid(p.pid, 0)

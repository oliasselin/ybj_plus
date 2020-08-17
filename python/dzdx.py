import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os
import subprocess
import sys

ts_plot=[100]#np.arange(0,152,1)#[66]#np.arange(0,181,10)#np.arange(141,201,1)#[50,100,150,200]

field = 1
contours=0
plot_slice=1
make_gif=0
show=1
quivers=1

field_title = ['dzdx']
field_name  = ['dz/dx']



#Run specification##############################################
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'N2_1e5_a'

constant_beta=0
ii=np.complex(1.j)
N2=1e-5
cor = 1.24e-4
f2=np.power(cor,2)
sigma_w=50.
zeta_max=2*np.pi/(100000*1.24e-4)  #Simply = k/f
grad_zeta=2.7e-9
U_scale = 0.5
u_0 = 0.1

hres=256
vres=256
timestep=0.1 #0.1 #Fraction of an inertial period between slices
################################################################



depth = 1000                   #Depth of interest (m)
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
kappa = k_over_sqrt2


Vprime = np.zeros((gp_depth,gp_xrange))          #V'(z,x')    --- along-stream velocity)
Zprime = np.zeros((gp_depth,gp_xrange))          #Zeta'(z,x') --- vertical vorticity (normalized by f) 
Zeta   = np.zeros((gp_depth,gp_xrange))          #Zeta'(z,x') --- vertical vorticity 
GradZp = np.zeros((gp_depth,gp_xrange))          #d Zeta'(z,x')/dx' --- gradient of vertical vorticity at y'=0                                    
ZZ     = np.zeros((gp_depth,gp_xrange))          #Depth vector (>0)
XP     = np.zeros((gp_depth,gp_xrange))          #x' vector (m)
for iz in range(gp_depth):
    Vprime[iz,:] = np.sqrt(2)*U_scale*np.cos(k_over_sqrt2*xp_m)
    Zprime[iz,:] = -zeta_max*np.sin(k_over_sqrt2*xp_m)
    GradZp[iz,:] = -grad_zeta*np.cos(k_over_sqrt2*xp_m)
    ZZ[iz,:]     = -dz*gpz[iz]
    XP[iz,:]     = xp_m[:]
    Zeta[iz,:]   = Zprime[iz,:]*cor

if(constant_beta==1):
    beta  = grad_zeta
else:
    beta  = GradZp


#############Plot parameters###############
colormap='RdBu_r' 
aspect=1.*gp_del/gp_depth
ncases=1
nts=1

nxticks=2
txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)
nyticks=5
tylabels=np.arange(0,depth+1,depth/nyticks)
ticksy_loc=np.arange(0,gp_depth,(gp_depth-1.)/nyticks)

wke_threshold = 1e-5
wke_levels=[1e-5, 1e-4, 1e-3, 1e-2]
lar_levels=[-0.02, -0.01, 0., 0.01, 0.02]
##########################################

Z, X = np.meshgrid(np.arange(0,gp_depth), np.arange(0,gp_xrange))

####### Tracing starting points ##########
xstart_list=[int(1*gp_del/4),int(gp_del/2),int(3*gp_del/4),gp_del,int(5*gp_del/4),int(3*gp_del/2),int(7*gp_del/4)]#[(int(gp_xrange/2)-20)]
zstart=1


#Create folder for plots if it doesn't exists                                                                                                                 
if not os.path.exists('plots/'+run+'/'+field_name[field-1]+'/'):
    os.makedirs('plots/'+run+'/'+field_name[field-1]+'/')

for ts in ts_plot:

    print "ts=",ts

    time = (2*np.pi/cor)*ts*timestep   #time is s
     
    
    dzdx = -np.power((3./2.)*cor*beta*ZZ/N2,(1./3.))
    
    
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
                im = ax.imshow(dzdx,cmap=colormap,aspect=aspect)

            if(quivers==1):

                #convert to give the right angle in the domain:
                stretch_x_m=xyrange_km*2*1000
                stretch_z_m=depth

                theta=dzdx#=np.arctan(dzdx)

                print(theta)

                cgx=stretch_z_m*np.cos(theta)
                cgz=stretch_x_m*np.sin(theta)
                
#                Q = ax.quiver(X[::8,::8],Z[::8,::8],cgx[::8,::8],cgz[::8,::8],color='k',pivot='mid', units='width')
                Q = ax.quiver(X[::8,::8],Z[::8,::8],cgx[::8,::8],cgz[::8,::8],color='k',pivot='mid', units='width')
                
                #cbar = ax.cax.colorbar(im)
            cbar = ax.cax.colorbar(im)#,ticks=lar_levels)# ticks=[1e-8, 1e-2])                                                                                           
            cbar = grid.cbar_axes[0].colorbar(im)
                
            if(contours==1):
                im=ax.contour(LAR,levels=lar_levels,colors='k')#,colors='k')

            time=ts*timestep
            time_title = '%.1f' % time
            ax.set_title(r''+field_title[field-1]+' (m/s)$^2$, $t =$ '+time_title+' inertial periods',fontsize=12)
            ax.text(-gp_del/5, gp_depth/2,r'Depth (m)',rotation='vertical',horizontalalignment='center',verticalalignment='center', fontsize=12)
            ax.text(gp_del, gp_depth+gp_depth/7,r"$x'$ (km)",rotation='horizontal',horizontalalignment='center',verticalalignment='center', fontsize=12)                

                
            if(show==1):
                plt.show()
            else:
                zeros_ts = (3-len(str(ts)))*'0'
                if(make_gif==1):
                    plt.savefig('plots/'+run+'/'+field_name[field-1]+'/'+field_name[field-1]+'_t'+zeros_ts+str(ts)+'.png',bbox_inches='tight')
                else:
                    plt.savefig('plots/'+run+'/'+field_name[field-1]+'/'+field_name[field-1]+'_t'+zeros_ts+str(ts)+'.eps',bbox_inches='tight')
                    
                    
if(make_gif==1):
    make_gif = 'convert -limit thread 1 -delay 1 -loop 0 plots/'+run+'/'+field_name[field-1]+'/*.png plots/'+run+'/'+field_name[field-1]+'/'+field_name[field-1]+'.gif'
    p = subprocess.Popen(make_gif, shell = True)
    os.waitpid(p.pid, 0)

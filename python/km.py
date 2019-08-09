import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import numpy as np
import os
import subprocess
import sys

ts_plot=np.arange(0,201,1)#[50,100,150,200]

show=0
plot_k=1
plot_m=0


#Run specification##############################################
scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'N2_1e5_dla'

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

Vprime = np.zeros((gp_depth,gp_xrange))          #V'(z,x')    --- along-stream velocity)
Zprime = np.zeros((gp_depth,gp_xrange))          #Zeta'(z,x') --- vertical vorticity 
for iz in range(gp_depth):
    Vprime[iz,:] = np.sqrt(2)*U_scale*np.cos(k_over_sqrt2*xp_m)
    Zprime[iz,:] = -zeta_max*np.sin(k_over_sqrt2*xp_m)

###########################################
#Quick and dirty plot uses one point
kmx0=int(4*gp_xrange/7)  #int(gp_xrange/2)#int(4*gp_xrange/7)
kmx1=int(5*gp_xrange/7)  #int(4*gp_xrange/7)#int(5*gp_xrange/7)
kmz0=int(1*gp_depth/5)
kmz1=int(2*gp_depth/5)

kmxa=[kmx0, int(kmx0+(kmx1-kmx0)/2), kmx1]
kmxa=np.concatenate((kmxa,kmxa,kmxa))

kmza=[kmz0, int(kmz0+(kmz1-kmz0)/2), kmz1]
kmza=np.repeat(kmza,3)

kmxa_km = kmxa*dx/1000-xyrange_km
kmza_m = -kmza*dz

xp0 = kmx0*dx/1000-xyrange_km
xp1 = kmx1*dx/1000-xyrange_km


gamma=0.5*grad_zeta*np.cos(k_over_sqrt2*((xp1+xp0)/2.)*1000)
beta =0.5*grad_zeta
ave_depth = (kmz0+kmz1)*dz/2

print "x' = [",xp0,",",xp1,"] km"
print "z' = [",-kmz0*dz,",",-kmz1*dz,"] m"
print ave_depth

kmt = np.zeros((len(ts_plot),3))       #Time series of t,k,m at point kmx,kmz
kmtp = np.zeros((len(ts_plot),3,len(kmxa)))       #Time series of t,k,m at point kmx,kmz
pred = np.zeros((len(ts_plot),5))             #Prediction for k: k = 0.5 (max_grad_zeta) t



#############Plot parameters###############
colormap='RdBu_r' 
aspect=0.6
ncases=1
nts=1

nxticks=2
txlabels=np.arange(-xyrange_km,xyrange_km+1,(2*xyrange_km)/nxticks)
ticksx_loc=np.arange(0,2*gp_del,(2*gp_del-1.)/nxticks)
nyticks=4
tylabels=np.arange(0,depth+1,depth/nyticks)
ticksy_loc=np.arange(0,gp_depth,(gp_depth-1.)/nyticks)

wke_threshold = 1e-5
wke_levels=[1e-5, 1e-4, 1e-3, 1e-2]
##########################################



####### Tracing starting points ##########
xstart_list=[int(1*gp_del/4),int(gp_del/2),int(3*gp_del/4),gp_del,int(5*gp_del/4),int(3*gp_del/2),int(7*gp_del/4)]#[(int(gp_xrange/2)-20)]
zstart=1


#Create folder for plots if it doesn't exists                                                                                                                 
if not os.path.exists('plots/'+run+'/km'):
    os.makedirs('plots/'+run+'/km')

for its,ts in enumerate(ts_plot):

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
        k = np.loadtxt(path_k)
        m = np.loadtxt(path_m)
        #Get the time series of k and m
        kmt[its,0] = ts*timestep  #Time in inertial periods
        kmt[its,1] = np.average(k[kmz0:kmz1,kmx0:kmx1])
        kmt[its,2] = np.average(m[kmz0:kmz1,kmx0:kmx1])
#        pred[its]   = 0.5*grad_zeta*(2*np.pi/cor)*ts*timestep
        pred[its,0] = ts*timestep  #Time in inertial periods
        pred[its,1]   = gamma*(2*np.pi/cor)*ts*timestep #0.5*grad_zeta*np.cos(k_over_sqrt2*(xp1-xp0)*1000)*(2*np.pi/cor)*ts*timestep
        pred[its,2]   = np.power((N2*gamma*gamma)/(3*cor*ave_depth),(1./3.))*(2*np.pi/cor)*ts*timestep
        pred[its,3]   =  beta*(2*np.pi/cor)*ts*timestep #0.5*grad_zeta*np.cos(k_over_sqrt2*(xp1-xp0)*1000)*(2*np.pi/cor)*ts*timestep
        pred[its,4]   = np.power((N2*beta*beta)/(3*cor*ave_depth),(1./3.))*(2*np.pi/cor)*ts*timestep


        for point in range(len(kmxa)):
            kmtp[its,0,point] = ts*timestep  #Time in inertial periods
            kmtp[its,1,point] = k[kmza[point],kmxa[point]]
            kmtp[its,2,point] = m[kmza[point],kmxa[point]]


    else:
        kmt = kmt[:its-1,:]
        kmtp=kmtp[:its-1,:,:]
        pred = pred[:its-1,:]




if(plot_k==1):
    plt.plot(kmt[:,0],kmt[:,1],color='k',label=r'$\bar{k}$',linewidth=2.)
    plt.plot(kmt[:,0],pred[:,3],'-b',label=r"$k = 0.5 \beta t$",linewidth=2.)
    plt.plot(kmt[:,0],pred[:,1],'-r',label=r"$k = - 0.5 \zeta_{x'}(\bar{x'}) t$",linewidth=2.)

    for point in range(len(kmxa)):
        plt.plot(kmtp[:,0,point],kmtp[:,1,point],color='gray',linewidth=0.3,label='_nolegend_')

elif(plot_m==1):
    plt.plot(kmt[:,0],kmt[:,2],color='k',label=r'$\bar{m}$',linewidth=3.)
    plt.plot(kmt[:,0],pred[:,4],'-b',label=r"$m = (N^2 \beta^2 / 12 f |\bar{z}|)^{1/3} t$",linewidth=3.)
    plt.plot(kmt[:,0],pred[:,2],'-r',label=r"$m = (N^2 \zeta_{x'}(\bar{x'})^2 / 12 f |\bar{z}|)^{1/3} t$",linewidth=3.)


    for point in range(len(kmxa)):
        plt.plot(kmtp[:,0,point],kmtp[:,2,point],color='gray',linewidth=0.3,label='_nolegend_')

#        plt.plot(kmtp[:,0,point],kmtp[:,1,point],label='x='+str(kmxa_km[point])+' km, z ='+str(kmza_m[point])+' m.')

plt.legend(loc='upper left',prop={'size': 9})
plt.xlabel(r'$t$ (inertial periods)')
plt.ylabel(r'Wavenumber ($m^{-1}$)')
    
if(show==1):
    plt.show()
else:
    zeros_ts = (3-len(str(ts)))*'0'
    if(plot_k==1):
        plt.savefig('plots/'+run+'/avek.eps',bbox_inches='tight')
    elif(plot_m==1):
        plt.savefig('plots/'+run+'/avem.eps',bbox_inches='tight')

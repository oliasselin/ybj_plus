import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

plot_phi=1
plot_u  =0


timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/double_gaussian/'
run = 'test_ml_sig33'

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run+'/'):
    os.makedirs('plots/'+run+'/')


#Defining the transect and converting distance in grid points
leif_dist = 9 #km 
leif_dist_x = leif_dist*np.cos(np.deg2rad(45)) #Projection on the x-axis assuming 45 deg angle
Dx = 100      #km domain size in the horizontal
gpdist = int(hres*leif_dist_x/Dx)
x0_km   = 0 #Distance away from the center
gpx0    = hres/2+int(hres*x0_km/Dx)

c_loc = gpx0-int(gpdist/2)          #location in x for the cyclone
a_loc = gpx0+gpdist-int(gpdist/2)   #location in x for the anti-cyclone


#Depth of the transect
v_ave = 1
leif_vave = 50 #m
Dz = 3000
gpave =  int(vres*leif_vave/Dz)
depth_top = 0 #in m              
gpdepth = int(vres*depth_top/Dz)

ticks = np.arange(0.,-12*np.pi-1.,-4*np.pi)
ticks_labels = ['$0$','$-4\pi$','$-8\pi$','$-12\pi$']


#Time range (time should be in fractions of inertial period. Fraction = timestep var above)
ts_min=0
ts_max=51#100

#Declare vars
ua = np.zeros(ts_max-ts_min)    
va = np.zeros(ts_max-ts_min)    
uc = np.zeros(ts_max-ts_min)    
vc = np.zeros(ts_max-ts_min)    
time = np.zeros(ts_max-ts_min)
stairs_c = np.zeros(ts_max-ts_min)
stairs_a = np.zeros(ts_max-ts_min)


#Suboptimal: load vertical slice
for ts in range(ts_min,ts_max):

    its=ts-ts_min

    spaces_ts = (3-len(str(ts)))*' '
    path_wke  = scratch_location+folder+run+'/output/slicev1'+spaces_ts+str(ts)+'.dat'
    path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
    path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'


    #Load the data file                                                                                                                  
    if os.path.isfile(path_wke):

        time[its]=ts*timestep #time in inertial periods

        print "Acquiring wave velocity data for t = ",time[its]," (index=",its,")"

        f_wke = np.loadtxt(path_wke)                 #Loads the full file as a 1-dim array                                                      
        f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                      
        f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                      

        g_wke   = np.rot90(np.reshape(f_wke,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
        g_lar   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
        g_lai   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                

        if(v_ave==1):
            lar_a = np.average(g_lar[gpdepth:(gpdepth+gpave),a_loc])          
            lai_a = np.average(g_lai[gpdepth:(gpdepth+gpave),a_loc])          
            
            lar_c = np.average(g_lar[gpdepth:(gpdepth+gpave),c_loc])          
            lai_c = np.average(g_lai[gpdepth:(gpdepth+gpave),c_loc])          
        else:
            lar_a = g_lar[gpdepth,a_loc]          
            lai_a = g_lai[gpdepth,a_loc]          
            
            lar_c = g_lar[gpdepth,c_loc]          
            lai_c = g_lai[gpdepth,c_loc]          

        ua[its] =  lar_a*np.cos(time[its]*2.*np.pi) + lai_a*np.sin(time[its]*2.*np.pi)
        va[its] = -lar_a*np.sin(time[its]*2.*np.pi) + lai_a*np.cos(time[its]*2.*np.pi)

        uc[its] =  lar_c*np.cos(time[its]*2.*np.pi) + lai_c*np.sin(time[its]*2.*np.pi)
        vc[its] = -lar_c*np.sin(time[its]*2.*np.pi) + lai_c*np.cos(time[its]*2.*np.pi)


        #Create a staircase function for both a and c points to make phi monotonuous
        if its >= 1:
            if np.arctan(vc[its]/uc[its]) > np.arctan(vc[its-1]/uc[its-1]):
                stairs_c[its]=stairs_c[its-1]+np.pi
            else:
                stairs_c[its]=stairs_c[its-1]

            if np.arctan(va[its]/ua[its]) > np.arctan(va[its-1]/ua[its-1]):
                stairs_a[its]=stairs_a[its-1]+np.pi
            else:
                stairs_a[its]=stairs_a[its-1]
        else:
            stairs_c[its]=0
            stairs_a[its]=0

#phi_final-phi_initial
dphia = np.arctan(va[-1]/ua[-1])-stairs_a[-1]  -  (np.arctan(va[0]/ua[0])-stairs_a[0])
dphic = np.arctan(vc[-1]/uc[-1])-stairs_c[-1]  -  (np.arctan(vc[0]/uc[0])-stairs_c[0])



zeta_eff_a = -2.*(1+dphia/(2.*np.pi*(time[-1]-time[0])))
zeta_eff_c = -2.*(1+dphic/(2.*np.pi*(time[-1]-time[0])))

f_eff_a = - dphia/(2.*np.pi*(time[-1]-time[0]))
f_eff_c = - dphic/(2.*np.pi*(time[-1]-time[0]))

#Last point:
print "A: (phi,time in IP)",dphia,time[-1]
print "C: (phi,time in IP)",dphic,time[-1]


#Just calculate vorticity at the a and c points
path_vort  = scratch_location+folder+run+'/output/slice2v7 0.dat'
f_vort = np.loadtxt(path_vort)
g_vort = np.rot90(np.reshape(f_vort,(vres,hres)),k=2) 

#vort_a = g_vort[0,a_loc]
#vort_c = g_vort[0,c_loc]

if(v_ave==1):
    vort_a = np.average(g_vort[gpdepth:(gpdepth+gpave),a_loc])
    vort_c = np.average(g_vort[gpdepth:(gpdepth+gpave),c_loc])

else:
    vort_a = g_vort[gpdepth,a_loc]
    vort_c = g_vort[gpdepth,c_loc]



print "Vorticity/f at the a and c spots:",vort_a,vort_c
print "Fit of phi  at the a and c spots:",zeta_eff_a,zeta_eff_c



labela = '%.2f' % f_eff_a
labelc = '%.2f' % f_eff_c

vorta = '%.2f' % float(vort_a/2.)
vortc = '%.2f' % float(vort_c/2.)

if(plot_phi==1):
    plt.plot(time,np.arctan(vc/uc)-stairs_c,'-*b',label=labelc+'$f$, $\zeta_c/2=$'+vortc+'$f$')
    plt.plot(time,np.arctan(va/ua)-stairs_a,'-*r',label=labela+'$f$, $\zeta_a/2=$'+vorta+'$f$')
    plt.xlabel(r'Time (inertial periods)')
    plt.ylabel(r'$\phi_u$ (rad)')
    plt.yticks(ticks,ticks_labels)
    plt.title(r'Angle of the wave velocity vector, $x_{cs}$ = '+str(x0_km)+' km, depth = '+str(depth_top)+' m')
    plt.ylim(ticks[-1],ticks[0])
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.legend(loc='best')
    plt.xlim(int(ts_min*timestep),int(ts_max*timestep))


#    plt.show()
    plt.savefig('plots/'+run+'/phi_d'+str(depth_top)+'.eps')

if(plot_u==1):
    plt.plot(time,uc,'-*b',label='u cyclonic')
    plt.plot(time,ua,'-*r',label='u anticyclonic')
    plt.xlabel(r'Time (inertial periods)')
    plt.ylabel(r'$u$ (m/s)')
    plt.title(r'Zonal velocity time series, $x_{cs}$ = '+str(x0_km)+' km, depth = '+str(depth_top)+' m' )
    plt.grid(color='k', linestyle='-', linewidth=0.25)
    plt.legend(loc='best')
    plt.xlim(int(ts_min*timestep),int(ts_max*timestep))
    
#    plt.show()
    plt.savefig('plots/'+run+'/u_d'+str(depth_top)+'.eps')

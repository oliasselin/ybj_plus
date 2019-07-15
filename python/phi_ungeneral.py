import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'test2_0.1out'

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'):
    os.makedirs('plots/')

leif_dist = 10 #km 
leif_dist_x = leif_dist*np.cos(np.deg2rad(45)) #Projection on the x-axis assuming 45 deg angle
Dx = 222      #km domain size in the horizontal
gpdist = int(hres*leif_dist_x/Dx)

c_loc = hres/2-int(gpdist/2)          #location in x for the cyclone
a_loc = hres/2+gpdist-int(gpdist/2)   #location in x for the anti-cyclone


ticks = np.arange(4*np.pi,-12*np.pi-1.,-4*np.pi)
ticks_labels = ['$4\pi$','$0$','$-4\pi$','$-8\pi$','$-12\pi$']

print "ticks and their labels: ",ticks,ticks_labels

ts=0
ts_max=51#100

ua = np.zeros(ts_max)    
va = np.zeros(ts_max)    
uc = np.zeros(ts_max)    
vc = np.zeros(ts_max)    
time = np.zeros(ts_max)
staircase = np.zeros(ts_max)


#Suboptimal: load vertical slice
for ts in range(0,ts_max):

    spaces_ts = (3-len(str(ts)))*' '
    path_wke  = scratch_location+folder+run+'/output/slicev1'+spaces_ts+str(ts)+'.dat'
    path_lar  = scratch_location+folder+run+'/output/slicev2'+spaces_ts+str(ts)+'.dat'
    path_lai  = scratch_location+folder+run+'/output/slicev3'+spaces_ts+str(ts)+'.dat'


    #Load the data file                                                                                                                  
    if os.path.isfile(path_wke):

        time[ts]=ts*timestep #time in inertial periods

        f_wke = np.loadtxt(path_wke)                 #Loads the full file as a 1-dim array                                                      
        f_lar = np.loadtxt(path_lar)                 #Loads the full file as a 1-dim array                                                      
        f_lai = np.loadtxt(path_lai)                 #Loads the full file as a 1-dim array                                                      

        g_wke   = np.rot90(np.reshape(f_wke,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
        g_lar   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
        g_lai   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                

        lar_a = g_lar[0,a_loc]          
        lai_a = g_lai[0,a_loc]          

        lar_c = g_lar[0,c_loc]          
        lai_c = g_lai[0,c_loc]          

        ua[ts] =  lar_a*np.cos(time[ts]*2.*np.pi) + lai_a*np.sin(time[ts]*2.*np.pi)
        va[ts] = -lar_a*np.sin(time[ts]*2.*np.pi) + lai_a*np.cos(time[ts]*2.*np.pi)

        uc[ts] =  lar_c*np.cos(time[ts]*2.*np.pi) + lai_c*np.sin(time[ts]*2.*np.pi)
        vc[ts] = -lar_c*np.sin(time[ts]*2.*np.pi) + lai_c*np.cos(time[ts]*2.*np.pi)

        staircase[ts] = np.pi*np.floor(2.*(time[ts]+0.25))    #to make the atan function smoother
        

#        wke[:,ts] = g_wke[0,:]          

        print 'ts=',ts
         
#plt.plot(time,uc,'-*b')
#plt.plot(time,ua,'-*r')



plt.plot(time,np.arctan(vc/uc)-staircase,'-*b')
plt.plot(time,np.arctan(va/ua)-staircase,'-*r')
plt.xlabel(r'Time (inertial periods)')
plt.ylabel(r'$\phi_u$ (rad)')
plt.yticks(ticks,ticks_labels)
plt.title('Angle of the wave velocity vector')
plt.ylim(ticks[-1],ticks[0])
plt.xlim(0,5)


#plt.show()
plt.savefig('plots/phi.png')

dphia = np.arctan(va[ts_max-1]/ua[ts_max-1])-staircase[ts_max-1]
dphic = np.arctan(vc[ts_max-1]/uc[ts_max-1])-staircase[ts_max-1]

zeta_eff_a = -2.*(1+dphia/(2.*np.pi*time[ts_max-1]))
zeta_eff_c = -2.*(1+dphic/(2.*np.pi*time[ts_max-1]))

#Last point:
print "A: (phi,time in IP)",np.arctan(va[ts_max-1]/ua[ts_max-1])-staircase[ts_max-1],time[ts_max-1]
print "C: (phi,time in IP)",np.arctan(vc[ts_max-1]/uc[ts_max-1])-staircase[ts_max-1],time[ts_max-1]


#Just calculate vorticity at the a and c points
path_vort  = scratch_location+folder+run+'_geo/output/slice2v7 0.dat'
f_vort = np.loadtxt(path_vort)
g_vort = np.rot90(np.reshape(f_vort,(vres,hres)),k=2) 

vort_a = g_vort[0,a_loc]
vort_c = g_vort[0,c_loc]

print "Vorticity/f at the a and c spots:",vort_a,vort_c
print "Fit of phi  at the a and c spots:",zeta_eff_a,zeta_eff_c







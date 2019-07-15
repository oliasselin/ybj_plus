import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

timestep=0.1 #0.1 #Fraction of an inertial period between slices

hres=256
vres=256

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'leif/'
run = 'attempt3_ro50'

plot_slice=0
colormap='RdBu_r' 
ncases=1
nts=1

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'+run+'/'):
    os.makedirs('plots/'+run+'/')

leif_dist = 9 #km 
leif_dist_x = leif_dist*np.cos(np.deg2rad(45)) #Projection on the x-axis assuming 45 deg angle
Dx = 100      #km domain size in the horizontal
gpdist = int(hres*leif_dist_x/Dx)

c_loc = hres/2-int(gpdist/2)          #location in x for the cyclone
a_loc = hres/2+gpdist-int(gpdist/2)   #location in x for the anti-cyclone


depth_top = 400 #in m
Dz = 3000      #ocean depth 
gpdepth = int(vres*depth_top/Dz)
leif_ave = 50  #top X m for average
gp_ave = 1#int(vres*leif_ave/Dz)

ticks = np.arange(4*np.pi,-12*np.pi-1.,-4*np.pi)
ticks_labels = ['$4\pi$','$0$','$-4\pi$','$-8\pi$','$-12\pi$']

print "ticks and their labels: ",ticks,ticks_labels

ts=0
ts_max=51#100

g_wke = np.zeros((vres,hres))    
g_lar = np.zeros((vres,hres))    
g_lai = np.zeros((vres,hres))    

ua = np.zeros(ts_max)    
va = np.zeros(ts_max)    
uc = np.zeros(ts_max)    
vc = np.zeros(ts_max)    
wke = np.zeros((hres,ts_max))    
time = np.zeros(ts_max)
staircase = np.zeros(ts_max)

#Just calculate vorticity at the a and c points
path_vort  = scratch_location+folder+run+'_geo/output/slice2v7 0.dat'
f_vort = np.loadtxt(path_vort)
g_vort = np.rot90(np.reshape(f_vort,(vres,hres)),k=2) 

vort_a = g_vort[0,a_loc]
vort_c = g_vort[0,c_loc]

print "Vorticity/f at the a and c spots:",vort_a,vort_c

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

        g_wke[:,:]   = np.rot90(np.reshape(f_wke,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
        g_lar[:,:]   = np.rot90(np.reshape(f_lar,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                
        g_lai[:,:]   = np.rot90(np.reshape(f_lai,(vres,hres)),k=2)        #Reshapes the array into a 2-d one  g(z,x,t)                                

        if(ts==0):
            print g_wke[:,0]

        lar_a = np.average(g_lar[gpdepth:(gpdepth+gp_ave),a_loc])          
        lai_a = np.average(g_lai[gpdepth:(gpdepth+gp_ave),a_loc])          

        lar_c = np.average(g_lar[gpdepth:(gpdepth+gp_ave),c_loc])          
        lai_c = np.average(g_lai[gpdepth:(gpdepth+gp_ave),c_loc])          

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
plt.title('Equivalent of F4.1-3 in cruise report, depth =',str(depth_top),' m')
plt.xticks(ticks,ticks_labels)

plt.show()

#Last point:
print "A: (phi,time in IP)",np.arctan(va[ts_max-1]/ua[ts_max-1])-staircase[ts_max-1],time[ts_max-1]
print "C: (phi,time in IP)",np.arctan(vc[ts_max-1]/uc[ts_max-1])-staircase[ts_max-1],time[ts_max-1]









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
        
        ax.get_xaxis().set_ticks([])
        ax.get_yaxis().set_ticks([])
        
        im = ax.imshow(g[:,:,0],cmap=colormap)
        
        cbar = ax.cax.colorbar(im)
        cbar = grid.cbar_axes[0].colorbar(im)
        
        plt.show()
        #    plt.savefig('plots/wke_m'+m+'.eps')


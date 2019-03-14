import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
import os

m_list=['8']
u_list=['0.25']
res_list=['_n256']
dt_list=['']
ts_list=['10','20','30','40','50']

colormap='RdBu_r'#colormap='coolwarm'#'seismic'
logscale=1
labels=['YBJ','BQ','YBJ$^+$']
ncases=len(labels)
nts=len(ts_list)
hres=64
Uw_scale=2.5e-5
WKE0 = 0.5*Uw_scale**2
max_allowed=100          #X times the initial WKE0 is allowed
enforce_max=1
max_enforced=100

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/'

#Create folder for plots if it doesn't exists
if not os.path.exists('plots/'):
    os.makedirs('plots/')

for m in m_list:
    for res in res_list:
        for u in u_list:
            for dt in dt_list:

                g = np.zeros((hres,hres,ncases*nts))
                bo_min=0.01
                bo_max=0.

                for its,ts in enumerate(ts_list):

                    run = 'm'+m+'_U'+u+res+dt

                    spaces_ts = (3-len(str(ts)))*' '
                    path_ybj  = scratch_location+folder+'YBJ/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
                    path_bo   = scratch_location+folder+'BO/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
                    path_ybjp = scratch_location+folder+'YBJp/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
                    
                    #Load the data file                                                                                                                  
                    if os.path.isfile(path_bo) and os.path.isfile(path_ybj) and os.path.isfile(path_ybjp):

                        f_ybj = np.loadtxt(path_ybj)                 #Loads the full file as a 1-dim array                                                      
                        f_bo = np.loadtxt(path_bo)                 #Loads the full file as a 1-dim array                                                              
                        f_ybjp = np.loadtxt(path_ybjp)                 #Loads the full file as a 1-dim array                                                    
                        
                        g[:,:,its*ncases]   = np.rot90(np.reshape(f_ybj,(hres,hres)) /WKE0)        #Reshapes the array into a 2-d one                                  
                        g[:,:,its*ncases+1] = np.rot90(np.reshape(f_bo,(hres,hres))  /WKE0)        #Reshapes the array into a 2-d one                                      
                        g[:,:,its*ncases+2] = np.rot90(np.reshape(f_ybjp,(hres,hres))/WKE0)        #Reshapes the array into a 2-d one 
                        
                        if(np.max(g[:,:,its*ncases+1])>bo_max):
                            bo_max = np.max(g[:,:,its*ncases+1])
                        

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

                    if(idx<nts):
                        if(idx==0):
                            ax.set_title('t = '+ts_list[idx],fontsize=12)    
                        elif(idx==(nts-1)):
                            ax.set_title(ts_list[idx]+' IP',fontsize=12)    
                        else:
                            ax.set_title(ts_list[idx],fontsize=12)    
                            
#                    if(idx<ncases):
#                        ax.set_title(labels[idx],fontsize=12)

                    ax.set_ylabel(labels[int(idx/nts)],fontsize=12,rotation=0,labelpad=20)

#                    if(np.mod(idx,ncases)):
#                        ax.set_title('IP = '+ts_list[int(idx/ncases)],loc='left')
                    if(logscale==1):
                        if(enforce_max==1):
                            im = ax.imshow(g[:,:,np.mod(idx,nts)*ncases+int(idx/nts)],cmap=colormap,norm=colors.LogNorm(vmin=bo_min,vmax=max_enforced) )
                        else:
                            im = ax.imshow(g[:,:,np.mod(idx,nts)*ncases+int(idx/nts)],cmap=colormap,norm=colors.LogNorm(vmin=bo_min,vmax=np.minimum(bo_max,max_allowed)) )
                    else:
                        im = ax.imshow(g[:,:,np.mod(idx,nts)*ncases+int(idx/nts)],cmap=colormap,vmin=bo_min,vmax=np.minimum(bo_max,max_allowed) )


                cbar = ax.cax.colorbar(im)
                if(logscale==1): #For some reason, no ticks appear in log unless I explicitly ask for them...
#                    cbar = grid.cbar_axes[0].colorbar(im, ticks=[0.01,0.1,1,10])
#                    cbar.ax.set_yticklabels(['10$^{-2}$','10$^{-1}$','10$^{0}$','10$^{1}$'])
                    cbar = grid.cbar_axes[0].colorbar(im, ticks=[0.01,0.1,1,10,100])
                    cbar.ax.set_yticklabels(['10$^{-2}$','10$^{-1}$','10$^{0}$','10$^{1}$','10$^{2}$'])
                else:
                    cbar = grid.cbar_axes[0].colorbar(im)

#                plt.title('WKE at the surface, run = '+run,y=1,fontsize=14)
#                plt.text(0,0.5,'20 IP')
#                plt.text(-0.1,0,'10 IP')
#                plt.text(-1,0.15,'50 IP')
#                fig.suptitle('WKE at the surface, m = '+m+', U = '+u,y=1,fontsize=16)
                plt.show()
                if(logscale==1):
                    plt.savefig('plots/m8_log.eps')
                else:
                    plt.savefig('plots/m8.eps')

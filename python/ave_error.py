import matplotlib.pyplot as plt
import numpy as np
import os

folder = 'YBJp/'
u_list = ['0.05','0.15','0.25']
ro_list= ['0.01','0.03','0.05']

factor=[8,2,1]

norm = 'L2'
folder = 'YBJp/'
xaxis_pad=0.3

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

fig = plt.figure(figsize=(6,3*len(u_list)))

for iu , u in enumerate(u_list):
    
    path_error = 'data/'+folder+'ave_error_'+norm+'_u'+u+'.dat'

    if os.path.isfile(path_error):

        error = np.loadtxt(path_error)

#        plt.subplots(nrows=len(u_list),ncols=1,num=iu,figsize=(6,3*len(u_list)))
        
        ax = fig.add_subplot(len(u_list),1,len(u_list)-iu)

        #plt.subplot(len(u_list),1,iu+1)
#        plt.plot(error[:,0],error[:,1]*100,'o-',color='b', mew='0', ms=10,label='$Ro$ = '+ro_list[iu]+', YBJ')
#        plt.plot(error[:,0],error[:,2]*100,'+-',color='g', mew='0', ms=10,label='$Ro$ = '+ro_list[iu]+', YBJ$^+$')
        plt.plot(error[:,0],error[:,1]*100,'-bo',linewidth=0.25, ms=7, label='$Ro$ = '+ro_list[iu]+', YBJ')
        plt.plot(error[:,0],error[:,2]*100,'-g+',linewidth=0.25, ms=10,label='$Ro$ = '+ro_list[iu]+', YBJ$^+$')
#        ax.grid(color='k', linestyle='-', linewidth=0.1)
        ax.yaxis.grid(color='k', linestyle='-', linewidth=0.1)
        plt.legend(loc='upper right',prop={'size': 10})        
        if(iu==0):
            plt.xlabel(r"$m'$")
        plt.ylabel(r'$\bar{E} \, (\%)$',rotation=90,labelpad=10)
        plt.xlim((1-xaxis_pad,12+xaxis_pad))
        plt.ylim((0,40/factor[iu]))
        plt.yticks(np.arange(0, 41./factor[iu], step=int(10./factor[iu])))
        plt.xticks(np.arange(1, 13, step=1))


    else:
        print "Couldn't find file ",path_error

#plt.title('Time-averaged error after 50 inertial periods')
plt.tight_layout()
#plt.savefig('plots/ave_error_'+norm+'_all.eps')
plt.show()

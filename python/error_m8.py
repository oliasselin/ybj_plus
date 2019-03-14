import matplotlib.pyplot as plt
import numpy as np
import os

m_list=['8']#['1','2','3','4','5','6','7','8','9','10','11','12']
u_list=['0.25']
res_list=['_n256']
dt_list=['']

folder = 'YBJp/'
norm = 'L2'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

tmax=50

for m in m_list:
    for res in res_list:
        for u in u_list:
            for dt in dt_list:

                run = 'm'+m+'_U'+u+res+dt
                path_error = 'data/'+folder+'error_'+norm+'_'+run+'.dat'
            
                if os.path.isfile(path_error):

                    print 'Creating a plot for run=',run

                    error = np.loadtxt(path_error)

                    fig, ax = plt.subplots(1,figsize=(6, 3))
                    plt.plot(error[:tmax+1,0],error[:tmax+1,1]*100,'-b',linewidth=0.7,markersize=3,label='YBJ')
                    plt.plot(error[:tmax+1,0],error[:tmax+1,2]*100,'-g',linewidth=0.7,markersize=5,label='YBJ$^+$')

#                    ax.grid(color='k', linestyle='-', linewidth=0.1)
                    ax.yaxis.grid(color='k', linestyle='-', linewidth=0.1)
                
                    plt.legend(loc='upper left', prop={'size': 10})        
                    plt.xlabel('Inertial periods')
#                    plt.ylabel('Error (%)',rotation=90,labelpad=10)
                    plt.ylabel(r'$E(t) \, (\%)$',rotation=90,labelpad=10)
                    plt.ylim((0,50))
                    plt.yticks(np.arange(0, 51., step=10))
                    plt.xlim((0,tmax))
                    plt.show()
                    plt.tight_layout()
#                    plt.savefig('plots/error_m'+m+'.eps')
                


import matplotlib.pyplot as plt
import numpy as np
import os

m_list=['1']#['1','2','3','4','5','6','7','8','9','10','11','12']
u_list=['0.25']
res_list=['_n256']
dt_list=['_moreslices']

folder = 'YBJp/'
norm='L2'
#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('plots/'):
    os.makedirs('plots/')

for m in m_list:
    for res in res_list:
        for u in u_list:
            for dt in dt_list:

                run = 'm'+m+'_U'+u+res+dt
                path_core = 'data/'+folder+'core_'+run+'.dat'            
                if os.path.isfile(path_core):

                    core = np.loadtxt(path_core)

                    fig, ax = plt.subplots(1,figsize=(6, 3))
#                    plt.plot(error[:,0],error[:,1]*100,'ob',markersize=3,label='YBJ')
#                    plt.plot(error[:,0],error[:,2]*100,'+g',markersize=5,label='YBJ$^+$')
                    plt.plot(core[:,0],core[:,1],'-b',linewidth=0.7)#,markersize=5,markevery=20,label='YBJ')
                    plt.plot(core[:,0],core[:,2],'-g',linewidth=0.7)#,markersize=7,markevery=20,label='YBJ$^+$')
                    plt.plot(core[:,0],core[:,3],'-r',linewidth=0.7)#,markersize=7,markevery=20,label='BQ')
#                    ax.yaxis.grid(color='k', linestyle='-', linewidth=0.1)                
#                    plt.legend(loc='upper right', prop={'size': 10})        
                    plt.xlabel('Inertial periods')
                    plt.ylabel(r'WKE$_c$',rotation=90,labelpad=10)
                    plt.ylim((0.8,2.2))
                    plt.yticks(np.arange(1, 2.1, step=0.2))
                    plt.xlim((0,50))
                    plt.tight_layout()
                    plt.show()



#                    plt.savefig('plots/core.eps')
                


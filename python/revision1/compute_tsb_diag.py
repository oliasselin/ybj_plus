#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np


m_list=['1']
u_list=['0.25']#['0.05','0.15','0.25']
res_list=['_n256']


hres=64
vres=256
Uw_scale=2.5e-5
U_scale =float(u_list[0])
#themax=int((hres/2)*(hres/2+0.5))   #bottom left quarter

more='' #'_moreslices'

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('data/'+folder+'tsb_uv/'):
    os.makedirs('data/'+folder+'tsb_uv/')


timestep = 1 #0.1  #Points per inertial period                                                                                                                                             
tmax = 501

L = 50000
f = 1e-4
Bu11=0.25

for m in m_list:
    for u in u_list:
        for res in res_list:

            Ro = float(u)/(f*L)
            eps = int(m)**2*Ro/Bu11
            freq= 2*np.pi*Bu11
            freqp=4*np.pi*Bu11/(2*int(m)**2+Bu11)
            print "Ro,eps,freq,freqp=",Ro,eps,freq,freqp

            WKE0 = 0.5*Uw_scale**2*(np.cos(int(m)*np.pi/vres)**2)

            run = 'm'+m+'_U'+u+res+more #'_moreslices'
            uv = np.zeros((tmax,3))
            uv_exp = np.zeros((tmax,3))

    
            for ts in range(tmax):

                spaces_ts = (3-len(str(ts)))*' '

                path_ybj_BR = scratch_location+folder+'YBJ/'+run+'/output/slicehtop2'+spaces_ts+str(ts)+'.dat'
                path_ybj_BI = scratch_location+folder+'YBJ/'+run+'/output/slicehtop3'+spaces_ts+str(ts)+'.dat'

                path_ybjp_BR = scratch_location+folder+'YBJp/'+run+'/output/slicehtop2'+spaces_ts+str(ts)+'.dat'
                path_ybjp_BI = scratch_location+folder+'YBJp/'+run+'/output/slicehtop3'+spaces_ts+str(ts)+'.dat'

                if os.path.isfile(path_ybj_BR) and os.path.isfile(path_ybj_BI) and os.path.isfile(path_ybjp_BR) and os.path.isfile(path_ybjp_BI):

                    print 'Calculating TSB diagnostic for the slice no ',ts,' for run = ',run


                    BR_ybj = np.loadtxt(path_ybj_BR)
                    BI_ybj = np.loadtxt(path_ybj_BI)

                    BR_ybjp = np.loadtxt(path_ybjp_BR)
                    BI_ybjp = np.loadtxt(path_ybjp_BI)

                    uv_ybj = np.mean(BR_ybj)**2 + np.mean(BI_ybj)**2
                    uv_ybjp = np.mean(BR_ybjp)**2 + np.mean(BI_ybjp)**2
                    

                    uv[ts,0]  = ts*timestep
                    uv[ts,1]  = uv_ybj/WKE0
                    uv[ts,2]  = uv_ybjp/WKE0

                    #WHAT IS EXPECTED FROM ANALYTIC SOLUTION???

                    #Expected strong-dispersion regime solution
#                    uv_exp[ts,0]  = ts*timestep
#                    uv_exp[ts,1]  = 1. + 2*eps*(1-np.cos(freq* wke_exp[ts,0]))    #YBJ
#                    uv_exp[ts,2]  = 1. + 2*eps*(1-np.cos(freqp*wke_exp[ts,0]))    #YBJp


                else:
                    break
                
            uv     = uv[:ts-1,:]
#            uv_exp = uv_exp[:ts-1,:]
            np.savetxt('data/'+folder+'uv_'+run+'.dat',uv)
#            np.savetxt('data/'+folder+'uv_exp_'+run+'.dat',uv_exp)

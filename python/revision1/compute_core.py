#!/usr/bin/env python                                                                                                                                                               
import os
import numpy as np


m_list=['1']
u_list=['0.25']#['0.05','0.15','0.25']
res_list=['_n256']


hres=64
vres=256
Uw_scale=2.5e-5
themax=int((hres/2)*(hres/2+0.5))   #bottom left quarter

more='_moreslices'

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('data/'+folder+'wke_blq/'):
    os.makedirs('data/'+folder+'wke_blq/')


timestep = 0.1  #Points per inertial period                                                                                                                                             
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

            run = 'm'+m+'_U'+u+res+more
            wke = np.zeros((tmax,4))
            wke_exp  = np.zeros((tmax,3))
            wke_exp2 = np.zeros((tmax,3))
    
            wke_bo_bq   = np.zeros((hres/2,hres/2))
            wke_ybj_bq  = np.zeros((hres/2,hres/2))
            wke_ybjp_bq = np.zeros((hres/2,hres/2))

            for ts in range(tmax):

                spaces_ts = (3-len(str(ts)))*' '
                path_bo   = scratch_location+folder+'BO/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
                path_ybj  = scratch_location+folder+'YBJ/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
                path_ybjp = scratch_location+folder+'YBJp/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'

                if os.path.isfile(path_bo) and os.path.isfile(path_ybj) and os.path.isfile(path_ybjp):
                    print 'Calculating WKE for the slice no ',ts,' for run = ',run

                    wke_bo   = np.loadtxt(path_bo)
                    wke_ybj  = np.loadtxt(path_ybj)
                    wke_ybjp = np.loadtxt(path_ybjp)
#                wke_bo   = np.rot90(np.reshape(np.loadtxt(path_bo),(hres,hres),order='F')/WKE0)
#                wke_ybj  = np.rot90(np.reshape(np.loadtxt(path_ybj),(hres,hres),order='F')/WKE0)
#                wke_ybjp = np.rot90(np.reshape(np.loadtxt(path_ybjp),(hres,hres),order='F')/WKE0) 

                    if(ts==0):
                        WKE0=wke_bo[0]
                        print "for m, u=",m,u,", WKE0=",WKE0,"expected=",0.5*(Uw_scale**2)*(np.cos(int(m)*np.pi/vres))**2


                    wke[ts,0]  = ts*timestep
                    wke[ts,1]  = wke_ybj[themax]/WKE0
                    wke[ts,2]  = wke_ybjp[themax]/WKE0
                    wke[ts,3]  = wke_bo[themax]/WKE0

                    #Expected strong-dispersion regime solution
                    wke_exp[ts,0]  = ts*timestep
                    wke_exp[ts,1]  = 1. + 2*eps*(1-np.cos(freq* wke_exp[ts,0]))    #YBJ
                    wke_exp[ts,2]  = 1. + 2*eps*(1-np.cos(freqp*wke_exp[ts,0]))    #YBJp

                    #Expected strong-dispersion regime solution, including the 1/2 lap Psi term of TSB's (3.34)
                    wke_exp2[ts,0]  = ts*timestep
                    wke_exp2[ts,1]  = 1. + (3./2.)*eps*(1-np.cos(freq* wke_exp2[ts,0]))    #YBJ
                    wke_exp2[ts,2]  = 1. + (3./2.)*eps*(1-np.cos(freqp*wke_exp2[ts,0]))    #YBJp

                else:
                    break
                
            wke      = wke[:ts-1,:]
            wke_exp  = wke_exp[:ts-1,:]
            wke_exp2 = wke_exp2[:ts-1,:]
            np.savetxt('data/'+folder+'core_'+run+'.dat',wke)
            np.savetxt('data/'+folder+'cexp_'+run+'.dat',wke_exp)
            np.savetxt('data/'+folder+'cexp2_'+run+'.dat',wke_exp2)

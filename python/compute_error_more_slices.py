#!/usr/bin/env python                                                                                                                                                                   
import os
import numpy as np

m_list=['1']#['1','2','3','4','5','6','7','8','9','10','11','12']
u_list=['0.25']
res='_n256'
dt='_moreslices'
bu=''#'_Bu0.083' #''

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp/'

#Create folder for plots if it doesn't exists                                                                                                                                           
if not os.path.exists('data/'+folder):
    os.makedirs('data/'+folder)


timestep = 0.1  #Points per inertial period                                                                                                                          
tmax = 501      #Number of IP to be included in plot




for u in u_list:

    ave_error_L1 = np.zeros((len(m_list),3)) #Average L1 error [m , error on YBJ, error on YBJ+] for a given u, as a function of m                       
    ave_error_L2 = np.zeros((len(m_list),3)) #Average L2 error [m , error on YBJ, error on YBJ+] for a given u, as a function of m

    for im,m in enumerate(m_list):
   
        run = 'm'+m+'_U'+u+res+dt+bu
        error_L1 = np.zeros((tmax,3))  #L1-Error time series for that run
        error_L2 = np.zeros((tmax,3))  #L2-Error time series for that run
        
        for ts in range(tmax):  #For every timestep, compute the L1- and L2-norm errors

            spaces_ts = (3-len(str(ts)))*' '
            path_bo = scratch_location+folder+'BO/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
            path_ybj = scratch_location+folder+'YBJ/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'
            path_ybjp= scratch_location+folder+'YBJp/'+run+'/output/slicehtop1'+spaces_ts+str(ts)+'.dat'

            if os.path.isfile(path_bo) and os.path.isfile(path_ybj) and os.path.isfile(path_ybjp):   #If those files exist, proceed to calculate
                print 'Calculating error for the slice no ',ts,' for run = ',run
                wke_bo   = np.loadtxt(path_bo)
                wke_ybj  = np.loadtxt(path_ybj)
                wke_ybjp = np.loadtxt(path_ybjp)

                #First columns of error is time in number of inertial periods
                error_L1[ts,0]  = ts*timestep
                error_L2[ts,0]  = ts*timestep

                #L1-norm error#  <|WKE_i - WKE_BO|>/<|WKE_i|>                                                                                     
                error_L1[ts,1]  = np.mean(np.absolute(wke_ybj  - wke_bo)) / np.mean(np.absolute(wke_bo))    #Same as (3.22) in TBS but with A -> wke                        
                error_L1[ts,2]  = np.mean(np.absolute(wke_ybjp - wke_bo)) / np.mean(np.absolute(wke_bo))    #Same as (3.22) in TBS but with A -> wke for YBJ+                  
 
                #L2-norm error#  sqrt{ <(WKE_i - WKE_BO)^2>/<WKE_BO^2> }                                                                                                        
                error_L2[ts,1]  = np.sqrt(np.mean(np.square(wke_ybj  - wke_bo)) / np.mean(np.square(wke_bo)))    #Same as (3.22) in TBS but with A -> wke              
                error_L2[ts,2]  = np.sqrt(np.mean(np.square(wke_ybjp - wke_bo)) / np.mean(np.square(wke_bo)))    #Same as (3.22) in TBS but with A -> wke for YBJ+       

            else:
                error_L1 = error_L1[:ts-1,:]
                error_L2 = error_L2[:ts-1,:]
                break

        if(ts > 1):
            np.savetxt('data/'+folder+'error_L1_'+run+'.dat',error_L1)
            np.savetxt('data/'+folder+'error_L2_'+run+'.dat',error_L2)
            
            ave_error_L1[im,0]=int(m)
            ave_error_L1[im,1]=np.average(error_L1[:,1])
            ave_error_L1[im,2]=np.average(error_L1[:,2])
            
            ave_error_L2[im,0]=int(m)
            ave_error_L2[im,1]=np.average(error_L2[:,1])
            ave_error_L2[im,2]=np.average(error_L2[:,2])

    #For each u, once you've run through all m's, print the average error for each norm
    np.savetxt('data/'+folder+'ave_error_L1_u'+u+bu+dt+'.dat',ave_error_L1)
    np.savetxt('data/'+folder+'ave_error_L2_u'+u+bu+dt+'.dat',ave_error_L2)


    print "E(YBJ)=",ave_error_L2[0,1],"E(YBJp)=",ave_error_L2[0,2]

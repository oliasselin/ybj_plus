import re
import os
import subprocess
import sys

m_list=['2']#['1','2','3','4','5','6','7','8','9','10','11','12']
u_list=['0.25']#['0.05','0.15','0.25']
eqn_list=['0','1']
eqn_label=['YBJ','YBJp']
eqn_label_to_change=['YBJp','YBJ']

scratch_location = '/oasis/scratch/comet/oasselin/temp_project/'
folder = 'YBJp'

if not os.path.exists(scratch_location+folder):
    os.makedirs(scratch_location+folder)

#log_file = open(scratch_location+folder+"/ybj_launch.log","w")
#sys.stdout = log_file

for eqn_num,eqn in enumerate(eqn_list):
    for m in m_list:
        for u in u_list:
        
            #print "Modifying m in parameters.f90..."
            change_m_parameter = "sed -i 's/mmm = [0-9]*/mmm = "+m+" /g' parameters.f90" 
            p = subprocess.Popen(change_m_parameter, shell = True)
            os.waitpid(p.pid, 0)

            #print "Modifying u in parameters.f90..."
            change_u_parameter = "sed -i 's/U_scale = [0-9]*.[0-9]*/U_scale = "+u+" /g' parameters.f90" 
            p = subprocess.Popen(change_u_parameter, shell = True)
            os.waitpid(p.pid, 0)

            #print "Modifying ybj_plus in parameters.f90..."
            change_ybj_parameter = "sed -i 's/ybj_plus = [0-9]*/ybj_plus = "+eqn+" /g' parameters.f90" 
            p = subprocess.Popen(change_ybj_parameter, shell = True)
            os.waitpid(p.pid, 0)

            #Read parameters.f90 from the original filee:                                                                                                                           
            with open ('parameters.f90', 'rt') as in_file:  # Open file for reading of text data.                                                    
                for line in in_file: # Store each line in a string variable "line"                                                                                                  
                    if line.find('mmm =') != -1:
                        mm=int(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
                        #print 'm =',mm
                        if(mm != int(m)): 
                            print "Oh no! m has no been changed!"
                        else:
                            print "Confirmation of parameter (1/3): m"
                    if line.find('ybj_plus') != -1:
                        ybjp=int(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
                        #print "ybj plus =",ybjp
                        if(ybjp != int(eqn)): 
                            print "Oh no! ybj_plus has no been changed!"
                        else:
                            print "Confirmation of parameter (2/3): y"
                    if line.find('U_scale =') != -1:
                        uu=float(re.findall(r"[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)[0])
                        print "uu =",uu
                        if(str(uu) != u): 
                            print "Oh no! u has no been changed!"
                        else:
                            print "Confirmation of parameter (3/3): u"
    
    
            #Modify the name of the run in the compiling script
            #print "Modifying the run's name"
            change_m_lcomet = "sed -i 's/m[0-9]*_/m"+m+"_/g' lcomet" 
            p = subprocess.Popen(change_m_lcomet, shell = True)
            os.waitpid(p.pid, 0)

            change_u_lcomet = "sed -i 's/U[0-9]*.[0-9]*_/U"+u+"_/g' lcomet"
            p = subprocess.Popen(change_u_lcomet, shell = True)
            os.waitpid(p.pid, 0)


            #print 'Equation = ',eqn," with no ",eqn_num,". Equation to change = ",eqn_label_to_change[eqn_num]
            string_to_change1='dir=/oasis/scratch/comet/oasselin/temp_project/'+folder+'/'+eqn_label_to_change[eqn_num]+'/'
            string_to_change2='dir=/oasis/scratch/comet/oasselin/temp_project/'+folder+'/'+eqn_label[eqn_num]+'/'
            change_dir = "sed -i 's!"+string_to_change1+"!"+string_to_change2+"!g' lcomet"
            p = subprocess.Popen(change_dir, shell = True)
            os.waitpid(p.pid, 0)


            #Sanity check: print the source and data directories
            with open ('lcomet', 'rt') as in_file:  # Open file for reading of text data.                                                                                            
                for line in in_file: # Store each line in a string variable "line"                                                                                                    
                    if line.find('sourdir=') != -1:
                        sourdir = str(line).rstrip("\n\r")
                    if line.find('datadir=') != -1:
                        datadir = str(line).rstrip("\n\r")

            #Compile lcomet
            print "Compiling and launching m = ",m,", u = ",u,", ybj_plus = ",eqn,"."
            print sourdir
            print datadir
            print '      '
            compile_lcomet = "./lcomet"
            p = subprocess.Popen(compile_lcomet, shell = True)
            os.waitpid(p.pid, 0)



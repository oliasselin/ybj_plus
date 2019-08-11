from make_a_gif import make_a_gif

scratch_location='/oasis/scratch/comet/oasselin/temp_project/'
home_location='/home/oasselin/'
folder = 'leif/'
run_list =['real_gaussian']# ['N2_1e5_a']
field_list = ['1']
#field_list = ['7']
sli_list = ['v']
#sli_list = ['htop','hmid','hbot']
##sli_list = ['htop','hmid','hbot']

print('Launching the gifmaker for '+str(len(field_list))+' field(s) in folder '+folder+' for '+str(len(run_list))+' run(s).')

for name in run_list:
    for field in field_list:
        for sli in sli_list:
            make_a_gif(folder=folder,run=name,sli=sli,field=field,nmax=200,fixed_cbrange='minmax',cbmin=-0.02,cbmax=0.02,hres=256,vres=256,timestep=0.1,scratch_location=scratch_location,home_location=home_location,U_scale=0.5,L_scale=100000/6.28,cor=1.24e-4)
#            make_a_gif(run=folder+name,sli=sli,field=field,nmax=100,fixed_cbrange='minmax',cbmin=-0.2,cbmax=0.2,hres=512,vres=512,timestep=0.01,scratch_location=scratch_location,home_location=home_location,U_scale=0.01,L_scale=222000/6.28)
#            make_a_gif(run=folder+name,sli=sli,field=field,nmax=100,fixed_cbrange='minmax',cbmin=-0.1,cbmax=0.1,hres=512,vres=512,timestep=0.01,scratch_location=scratch_location,home_location=home_location,U_scale=0.01,L_scale=222000/6.28)

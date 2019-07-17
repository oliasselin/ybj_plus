import os
import subprocess
import sys



run = 'escape_aug14/w_L0_G1'
sli = 'v'   #htop,hmid,hbot,v
field = '2'

if field=='1':
    field_name = 'WKE'
if field=='2':
    field_name = 'Re-LA'
if field=='3':
    field_name = 'Im-LA'
if field=='4':
    field_name = 'WPE'

hres = 128
vres = 128

timestep=0.01

U_scale = 0.1
L_scale = 400000/(3.14159*2.)

tau_e = L_scale/(U_scale)/3600 #eddy turnover time in hours

delay = 1 #In cs

fixed_cbrange='minmax'     #0: free, 1: set min only, 2: set max only, 3: set both max and min 
cbmin = -0.005
cbmax = 0.005#1e-5

nmax = 350 #Maximum number of slices


for k in range(0,nmax):
    
    if k<10: 
        path_file_xy  = '/scratch/05518/oasselin/'+run+'/output/slice'+sli+field+'  '+str(k)+'.dat'
    if (k<100 and k>9):
        path_file_xy  = '/scratch/05518/oasselin/'+run+'/output/slice'+sli+field+' '+str(k)+'.dat'
    if k>=100:
        path_file_xy  = '/scratch/05518/oasselin/'+run+'/output/slice'+sli+field+str(k)+'.dat'

    if os.path.isfile(path_file_xy): 

        if k==0:
            png_dir='/scratch/05518/oasselin/'+run+'/temp/'
            if not os.path.exists(png_dir):
                os.makedirs(png_dir)

        if k<10:
            output_file = png_dir+'slice00'+str(k)+'.png'
        if (k<100 and k>9):
            output_file = png_dir+'slice0'+str(k)+'.png'
        if k>=100:
            output_file = png_dir+'slice'+str(k)+'.png'
        
        time = "{0:.2f}".format(timestep*k)
        time_hours = "{0:.2f}".format(timestep*k*tau_e)
        xlabel = 't='+str(time)+' tau_e ('+time_hours+' hours)'
        tt = field_name+' '+sli+' '+run


        if fixed_cbrange == 'minmax':
            cbrange='['+str(cbmin)+':'+str(cbmax)+']'
        elif fixed_cbrange == 'max':
            cbrange='[*:'+str(cbmax)+']'
        elif fixed_cbrange == 'min':
            cbrange='['+str(cbmin)+':*]'
        else:
            cbrange='[*:*]'


        xrange = '[0:'+str(hres-1)+']'
        if sli =='v':
            yrange = '[0:'+str(vres-1)+']'
        else:
            yrange = '[0:'+str(hres-1)+']'




        gnuplot_command = "gnuplot -e \"set output '"+output_file+"'; set xrange "+xrange+"; set yrange "+yrange+"; set cbrange "+cbrange+"; set title '"+tt+"'; set xlabel '"+xlabel+"'; filename = '"+path_file_xy+"'\" slice.gnu"

        p = subprocess.Popen(gnuplot_command, shell = True)
        os.waitpid(p.pid, 0)

gif_dir='/home1/05518/oasselin/gif/'+run

if not os.path.exists(gif_dir):
    os.makedirs(gif_dir)

make_gif = 'convert -delay '+str(delay)+' -loop 0 '+png_dir+'*.png '+gif_dir+'/'+field_name+'_'+sli+'.gif'

p = subprocess.Popen(make_gif, shell = True)
os.waitpid(p.pid, 0)


delete_content = 'rm '+png_dir+'*.png'
p = subprocess.Popen(delete_content, shell = True)
os.waitpid(p.pid, 0)

delete_folder = 'rmdir '+png_dir
p = subprocess.Popen(delete_folder, shell = True)
os.waitpid(p.pid, 0)

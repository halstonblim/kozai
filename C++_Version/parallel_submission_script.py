import multiprocessing
import subprocess
import numpy as np
import os
import argparse

flags = ["--quad","--oct","--peri","--rad","--ignore_gsl","--cross_lim","--outerperi","--print_gsl"]

## The filenames for the triples we'll be integrating

parser = argparse.ArgumentParser(description="Submit parallel jobs kozai.sh")
parser.add_argument("--cross_lim",action='store_true',default=False,help="turn interaction flag on")
parser.add_argument("--output_folder",type=str)
parser.add_argument("--input_file",type=str)
args = parser.parse_args()

param = ["--m1","--m2","--m3","--a1","--a2","--e1","--e2","--g1","--g2",
		 "--omega1","--omega2","--inc","--time","--dt"]

## Where the actual kozai program is
prog = "/dir/here"
output_folder = args.output_folder

## Some constants
age_univ = 1.38e10
int_time = age_univ
int_step = age_univ / 10**6

## Flags for integration
quad = True 
octp = True 
peri = True
outerperi = True
radr = True
igsl = True
pgsl = True
iinteract = args.cross_lim

## Set these flags into command line arguments
flags_str = []
if quad:
    flags_str.append(flags[0])
if octp:
    flags_str.append(flags[1])
if peri:
    flags_str.append(flags[2])
if outerperi:
    flags_str.append(flags[6])
if radr:
    flags_str.append(flags[3])
if igsl:
    flags_str.append(flags[4])
if iinteract:
    flags_str.append(flags[5])
if pgsl:
    flags_str.append(flags[7])

Tt = param[12]+' '+str(int_time)
dt = param[13]+' '+str(int_step)
def integrate_with_kozai_code(triple_str):
	## Must be check_call instead of Popen
	## the former just runs the command on the existing thread,
	## whereas the later actually opens a new kernel process

	## If you do that, you'll keep opening processes with the multiprocessing
	## pool until the machine kernel breaks...
  try:
    kozai_str = subprocess.run(triple_str,shell=True,check=True)
    fout  = open(output_folder + '/startedruns.txt','a')
    fout.write(triple_str + '\n')
    fout.close()
  except subprocess.CalledProcessError as e:
    fout  = open(output_folder + '/failedruns.txt','a')
    fout.write(triple_str + '\n')
    fout.close()

## Initialize the multiprocess pool; note this must come after the above
## definition, so each of the processes knows about it
p = multiprocessing.Pool()

all_run_strings = []

## load the file with the triples and make a new directory for the output
filename = args.input_file

if not os.path.exists(args.output_folder):
  os.makedirs(args.output_folder)

basename = args.output_folder + '/' + args.input_file.rstrip('.dat')
data = np.loadtxt(filename,ndmin=2)
print(np.shape(data))
tI = 1
for d in data:
  tI += 1
  m1 = param[0]+' '+str(d[0]); m2 = param[1]+' '+str(d[1])
  m3 = param[2]+' '+str(d[2])
  a1 = param[3]+' '+str(d[3]); a2 = param[4]+' '+str(d[4])
  e1 = param[5]+' '+str(d[5]); e2 = param[6]+' '+str(d[6])
  g1 = param[7]+' '+str(d[7]); g2 = param[8]+' '+str(d[8])
  if d[9] < 180:
    om2 = d[9] + 180
  else:
    om2 = d[9] - 180
  omega1 = param[9] +' '+str(d[9]); omega2 = param[10]+' '+str(om2)
  inc = param[11]+' '+str(d[11])

  err_name=basename+'-out.'+str(tI); 
  out_name=basename+'-inf.'+str(tI)

  triple_tI_base = [prog,Tt,dt,m1,m2,m3,a1,a2,e1,e2,g1,g2,omega1,omega2,inc]
  triple_tI = list(triple_tI_base); triple_tI.extend(flags_str)

  triple_str = ''
  for param_i in triple_tI:
    triple_str = triple_str+' '+param_i
  triple_str = triple_str[1:]
  triple_str += ' > '+out_name+' 2> '+err_name
  if(os.path.isfile(out_name) == False):
    all_run_strings.append(triple_str)

## Finally, map the calls to the kozai code out to the process pool
if __name__ == "__main__":
  print(all_run_strings[0])
  p.map(integrate_with_kozai_code,all_run_strings)

import sys
import os
import argparse
from math import *

with open("command.txt", "w") as of:
    of.write(" ".join(["python"]+sys.argv))

parser = argparse.ArgumentParser()

parser.add_argument("-d", "--inputdir", type=str, help="InputDir", required=True)
parser.add_argument("-c", "--condadir", type=str, help="CondaDir", required=True)
parser.add_argument("-e", "--condaenv", type=str, help="CondaEnv", required=True)
parser.add_argument("-f", "--config", type=str, help="Config", required=True)
parser.add_argument("-i", "--numiter", type=str, help="NumIter", required=True)
parser.add_argument("-n", "--numbound", type=str, help="NumBound", required=True)
parser.add_argument("-s", "--scalesig", type=str, help="ScaleSig", required=False)
parser.add_argument("-b", "--scalebkg", type=str, help="ScaleBkg", required=False)
parser.add_argument("-l", "--log", type=str, help="Log", required=True)
parser.add_argument("-B", "--bound", type=str, help="Boundaries", required=False)
parser.add_argument("-q", "--queue", type=str, help="Condor queue", default="longlunch", required=True)
args = parser.parse_args()

if not os.path.isdir('error'): os.mkdir('error') 
if not os.path.isdir('output'): os.mkdir('output') 
if not os.path.isdir('log'): os.mkdir('log') 


# Prepare condor jobs
condor = '''executable              = run_script.sh
output                  = output/strips.$(ClusterId).$(ProcId).out
error                   = error/strips.$(ClusterId).$(ProcId).err
log                     = log/strips.$(ClusterId).log
transfer_input_files    = run_script.sh
transfer_output_files   = ""
on_exit_remove          = (ExitBySignal == False) && (ExitCode == 0)
periodic_release        = (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > (60*60))

+AccountingGroup        = "group_u_CMS.CAF.ALCA" 
+JobFlavour             = "{queue}"
queue arguments from arguments.txt

'''

condor = condor.replace("{queue}", args.queue)

command = './optimize.py -i ./{CONFIG} -o {OUTNAME} --num-iter={ITER} --num-bounds={NUM} --log {LOG}'
if args.scalesig:
   command = command + ' --lumi-scale {SCALESIG}'
if args.scalebkg:
   command = command + ' --lumi-scale-bkg {SCALEBKG}'   
if args.bound:
   command = command + ' -b {BOUND}'   
#print "Command:",command

script = '''#!/bin/sh -e
eval "$(conda shell.bash hook)"
cd {INDIR}
conda config --prepend envs_dirs {CONDADIR}/envs
conda config --prepend pkgs_dirs {CONDADIR}/pkgs
conda activate {CONDADIR}/envs/{CONDAENV}
{COMMAND}

echo -e "DONE";
'''

log_output = args.log
outname = log_output.split('/')[-1]
outname = outname.replace('logs_','')
outname = outname.replace('.txt','')

script = script.replace("{COMMAND}", command)
script = script.replace("{CONDADIR}", args.condadir)
script = script.replace("{CONDAENV}", args.condaenv)
script = script.replace("{INDIR}", args.inputdir)
script = script.replace("{CONFIG}", args.config)
script = script.replace("{OUTNAME}", outname)
script = script.replace("{ITER}", args.numiter)
script = script.replace("{NUM}", args.numbound)
script = script.replace("{LOG}", args.log)
if args.scalesig:
   script = script.replace("{SCALESIG}", args.scalesig)
if args.scalebkg:
   script = script.replace("{SCALEBKG}", args.scalebkg)
if args.bound:
   script = script.replace("{BOUND}", args.bound)
    
with open("condor_job.txt", "w") as cnd_out:
    cnd_out.write(condor)

with open("run_script.sh", "w") as rs:
    rs.write(script)

#os.system("condor_submit condor_job.txt")





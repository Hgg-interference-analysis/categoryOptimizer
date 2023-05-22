import sys
import os
import argparse
import random
import commands
from math import *

with open("command.txt", "w") as of:
    of.write(" ".join(["python"]+sys.argv))

parser = argparse.ArgumentParser()

parser.add_argument("-o", "--outputdir", type=str, help="Outputdir", required=True)
parser.add_argument("-c", "--cmssw", type=str, help="CMSSW", required=True)
parser.add_argument("-q", "--queue", type=str, help="Condor queue", default="longlunch", required=True)
args = parser.parse_args()

CWD = os.getcwd()
OUTPUT = args.outputdir
if not os.path.isdir('error'): os.mkdir('error') 
if not os.path.isdir('output'): os.mkdir('output') 
if not os.path.isdir('log'): os.mkdir('log') 
if not os.path.isdir(OUTPUT+'/.conda'): os.mkdir(OUTPUT+'/.conda') 
if not os.path.isdir(OUTPUT+'/.conda/envs'): os.mkdir(OUTPUT+'/.conda/envs') 
if not os.path.isdir(OUTPUT+'/.conda/pkgs'): os.mkdir(OUTPUT+'/.conda/pkgs') 


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


script = '''#!/bin/sh -e
cd {CMSSW}/src/
echo -e "evaluate"
eval `scramv1 ru -sh`

cd {CWD}/
cd ../
chmod 755 -R {OUTPUT}/.conda
conda config --prepend envs_dirs {OUTPUT}/.conda/envs
conda config --prepend pkgs_dirs {OUTPUT}/.conda/pkgs
conda env create -f environment.yml

echo -e "DONE";
'''

script = script.replace("{CMSSW}", args.cmssw)
script = script.replace("{CWD}", CWD)
script = script.replace("{OUTPUT}", args.outputdir)
    
with open("condor_job.txt", "w") as cnd_out:
    cnd_out.write(condor)

with open("run_script.sh", "w") as rs:
    rs.write(script)

#os.system("condor_submit condor_job.txt")





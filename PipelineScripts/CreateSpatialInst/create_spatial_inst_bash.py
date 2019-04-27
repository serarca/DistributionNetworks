import pickle
import sys
import os
import collections
import json


folder_project = sys.argv[1]
results_folder_ts = "Results/" + sys.argv[2]
dist_type = sys.argv[3]
delta = sys.argv[4]
penalty = sys.argv[5]

if penalty != "no_penalty":
    penalty_factor = sys.argv[6]
else:
    penalty_factor = 0.0

gamma = 20000
time_dual = "9:00:00"
time_primal = "6:00:00"
memory_dual = "8G"
memory_primal = "16G"



folder_path = results_folder_ts+"bash_c/spatial/"
submission = results_folder_ts+"bash_c/bash_statusquo.sh"
cluster_path = "~/vrp_project/"
local_path = folder_project

os.chdir(folder_project)
sys.path.insert(0, 'Code/VRPEngine/pyCode')
sys.path.insert(0, 'Code/VRPEngine/C++Engine')
sys.path.insert(0, 'Code/VRPEngine/pyCode/tsp')

import VRPClass
import NetworkClass
infile = open(results_folder_ts+"GVRP/GVRP.p",'rb')
GVRP = pickle.load(infile)
infile.close()
GVRP.dist_type = dist_type

infile = open(results_folder_ts+"GVRP/Network.p",'rb')
Network = pickle.load(infile)
infile.close()


# Create bash scripts for the cluster
header_spatial_dual = """#!/bin/bash
#
#SBATCH --job-name=test
#
#SBATCH --time=%s
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

module load gurobi/8.0.1_py27
cd ~/vrp_project/Code/VRPEngine/C++Engine
srun ./main ~/vrp_project/%s spatial/spatial_day_%s %s -p normal,dev,gpu
"""

# Create bash scripts for the cluster
header_spatial_primal = """#!/bin/bash
#
#SBATCH --job-name=test
#
#SBATCH --time=%s
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16G

module load gurobi/8.0.1_py27
cd ~/vrp_project/Code/VRPEngine/C++Engine
srun ./primal ~/vrp_project/%s spatial/spatial_day_%s %s %s %s -p normal,dev,gpu\n
"""


folder_path = results_folder_ts+"bash_c/spatial_dual/"
submission = results_folder_ts+"bash_c/bash_spatial_dual.sh"
cluster_path = "~/vrp_project/"

open(submission, 'w').close()
with open(submission, 'a') as file:
    file.write('cd ~/vrp_project/%sbash_c\n'%results_folder_ts)


for day in GVRP.days_list:
    instance = "spatial_day_"+str(day)
    path = folder_path + instance +".sh"
    open(path, 'w').close()
    with open(path, 'a') as file:
        file.write(header_spatial_dual % (time_dual,results_folder_ts,str(day),str(penalty_factor)))
    with open(submission, 'a') as file:
        file.write('sbatch ' + cluster_path + results_folder_ts + 'bash_c/'+'spatial_dual/'+instance + '.sh' + '\n')
    

# Spatial primal
folder_path = results_folder_ts+"bash_c/spatial_primal/"
submission = results_folder_ts+"bash_c/bash_spatial_primal_"+delta+".sh"
cluster_path = "~/vrp_project/"

open(submission, 'w').close()
with open(submission, 'a') as file:
    file.write('cd ~/vrp_project/%sbash_c\n'%results_folder_ts)

for day in GVRP.days_list:
    instance = "spatial_day_"+str(day)+ "_" +str(delta)
    path = folder_path + instance +".sh"
    open(path, 'w').close()
    with open(path, 'a') as file:
        file.write(header_spatial_primal % (time_primal,results_folder_ts,str(day),str(delta), str(gamma), str(penalty_factor)))
    with open(submission, 'a') as file:
        file.write('sbatch ' + cluster_path + results_folder_ts + 'bash_c/'+'spatial_primal/'+instance + '.sh' + '\n')


# Spatial VRP instances
submission = results_folder_ts+"bash_c/bash_spatial_instances.sh"

open(submission, 'w').close()
for day in GVRP.days_list:
    with open(submission, 'a') as file:
        file.write('python2 %s/Code/PipelineScripts/CreateSpatialInst/create_spatial_inst.py %s %s %s %s %s'%(local_path,local_path, sys.argv[2], dist_type, penalty, str(day)) + '\n')






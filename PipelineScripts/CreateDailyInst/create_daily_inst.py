import pickle
import sys
import os
import collections
import json


folder_project = sys.argv[1]
results_folder_ts = "Results/" + sys.argv[2]
dist_type = sys.argv[3]


folder_path = results_folder_ts+"bash_c/statusquo/"
submission = results_folder_ts+"bash_c/bash_statusquo.sh"
cluster_path = "~/vrp_project/"

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
header_sq = """#!/bin/bash
#
#SBATCH --job-name=test
#
#SBATCH --time=%s
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

cd ~/vrp_project/Code/VRPEngine/C++Engine
module load gurobi/8.0.1_py27

"""
line_sq_dual = "srun ./main ~/vrp_project/%s statusquo/statusquo_cluster_%s_day_%s 0.0 -p normal,dev,gpu\n"
line_sq_primal = "srun ./primal ~/vrp_project/%s statusquo/statusquo_cluster_%s_day_%s %s %s 0.0 -p normal,dev,gpu\n"

open(submission, 'w').close()
with open(submission, 'a') as file:
    file.write('cd ~/vrp_project/%sbash_c\n'%results_folder_ts)

for cluster in GVRP.middlemen.keys():
    instance = "statusquo_cluster_"+str(cluster)
    path = folder_path + instance +".sh"
    open(path, 'w').close()
    with open(path, 'a') as file:
        file.write(header_sq % ("1:00:00"))
    for day in GVRP.days_list:
        fs = [f for f in GVRP.farmers.values() if (day in f.pickups and f.cluster_id_sq == cluster)]
        if (len(fs)>0):
            with open(path, 'a') as file:
                file.write(line_sq_dual % (results_folder_ts, str(cluster),str(day)))
                file.write(line_sq_primal % (results_folder_ts, str(cluster),str(day), "2000", "20000"))
    with open(submission, 'a') as file:
        file.write('sbatch ' + cluster_path + results_folder_ts + 'bash_c/'+'statusquo/'+instance + '.sh' + '\n')
    




# Create daily instances
n_farmers = collections.defaultdict(lambda : 0)

for cluster in GVRP.middlemen.keys():
    for day in GVRP.days_list:
        # Take all mills that belong to that middleman
        ms = [m for m in GVRP.mills.values() if (m.node_id in GVRP.middlemen[cluster].mill_codes_sq)]
        # Take all farmers with corresponding pickup day and cluster
        fs = [f for f in GVRP.farmers.values() if (day in f.pickups and f.cluster_id_sq == cluster)]
        # Take all trucks of a given cluster
        ts = [t for t in GVRP.trucks.values() if (t.cluster_id == cluster)]
        metadata = {
            "type":"status-quo",
            "day": day,
            "cluster": cluster
        }
        if (len(fs)>0):
            n_farmers[cluster]+=len(fs)
            vrp = VRPClass.VRP.create_VRP_instance(ms, fs, ts, metadata, dist_type=dist_type, network = Network)
            assert(vrp.feasibility())
            json_dict = VRPClass.c_json_parser(vrp)
            with open( results_folder_ts+"instances/statusquo/statusquo_"+"cluster_"+str(cluster)+"_day_"+str(day)+".json", "wb" ) as fp:
                json.dump(json_dict, fp, indent=4)
    print("Cluster %d: %d farmers"%(cluster,n_farmers[cluster]))

import pickle
import sys
import os
import collections
import json
from random import sample 


folder_project = sys.argv[1]
results_folder_ts = "Results/" + sys.argv[2]
dist_type = "road_time"
penalty = "village"

day = 0



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



# Create spatial instances
n_farmers = collections.defaultdict(lambda : 0)

# Take all mills
ms = [m for m in GVRP.mills.values()]
# Take all farmers with corresponding pickup day and cluster
fs = [f for f in GVRP.farmers.values() if (day in f.pickups)]
# Take all trucks
ts = [t for t in GVRP.trucks.values()]
metadata = {
    "type":"spatial",
    "day": day
}
if (len(fs)>0):
    #Debug
    if True:
        fs = [fs[i] for i in sample(list(range(len(fs))),3)]
        ts = [ts[i] for i in sample(list(range(len(ts))),3)]
        ms = [GVRP.mills['SKIP']]
    n_farmers[day]+=len(fs)
    print("Constructing VRP")
    penalty_type = None
    if penalty != "no_penalty":
        penalty_type = penalty
    vrp = VRPClass.VRP.create_VRP_instance(ms, fs, ts, metadata, dist_type=dist_type, network = Network, penalty_type = penalty_type)
    json_dict = VRPClass.c_json_parser(vrp)
    with open( results_folder_ts+"instances/spatial/tmp_spatial_"+"day_"+str(day)+".json", "wb" ) as fp:
        json.dump(json_dict, fp, indent=4)
print("Day %d: %d farmers"%(day,n_farmers[day]))

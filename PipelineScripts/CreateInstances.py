
# coding: utf-8

# In[1]:
from __future__ import print_function

import datetime
from dateutil import parser
import math
from dateutil.parser import parse
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
from shapely.geometry import Point
import json
import ast
import sys
import pickle
import os
import time
import collections

folder_project = "/Users/sergiocamelo/Dropbox/Sergio-Joann/"
os.chdir(folder_project)

sys.path.insert(0, 'Code/VRPEngine/pyCode')
sys.path.insert(0, 'Code/VRPEngine/C++Engine')
sys.path.insert(0, 'Code/VRPEngine/pyCode/tsp')

import solver as solver
import distances as distances
import VRPClass


# In[2]:

days_analysis = 14
ts = sys.argv[1]
data_folder = 'StandardizedData/'
results_folder_ts = 'Results/'+ ts + "/" 
# The folder which will be used to write the bash scripts


# In[3]:

# Add metadata file
metadata = {
    "time":ts
}
with open( results_folder_ts+"json_metadata.json", "wb" ) as fp:
    json.dump(metadata, fp)


# In[4]:

# Create results folders
try:
    os.mkdir(results_folder_ts+'tables_for_report/')
    os.mkdir(results_folder_ts+'data_cleaning_results/')
    os.mkdir(results_folder_ts+'vrps/')
    os.mkdir(results_folder_ts+'vrps/spatial/')
    os.mkdir(results_folder_ts+'vrps/temporal/')
    os.mkdir(results_folder_ts+'vrps/daily/')
    os.mkdir(results_folder_ts+'instances/')
    os.mkdir(results_folder_ts+'instances/spatial/')
    os.mkdir(results_folder_ts+'instances/temporal/')
    os.mkdir(results_folder_ts+'instances/daily/')
    os.mkdir(results_folder_ts+'low_cost_routes/')
    os.mkdir(results_folder_ts+'low_cost_routes/spatial/')
    os.mkdir(results_folder_ts+'low_cost_routes/temporal/')
    os.mkdir(results_folder_ts+'low_cost_routes/daily/')
    os.mkdir(results_folder_ts+'solution_routes/')
    os.mkdir(results_folder_ts+'solution_routes/spatial/')
    os.mkdir(results_folder_ts+'solution_routes/temporal/')
    os.mkdir(results_folder_ts+'solution_routes/daily/')
    os.mkdir(results_folder_ts+'bash/')
    os.mkdir(results_folder_ts+'bash/cluster/')
except:
    print("files created already")


# In[5]:

open(results_folder_ts + '/errors.txt', 'a').close()
f = open(results_folder_ts + '/errors.txt', 'w+')


# In[6]:

# Load data
dim_pickups = pd.read_csv(results_folder_ts+"dim_pickups.csv")
dim_pickups['pickups'] = dim_pickups['pickups'].map(eval)
dim_middlemen = pd.read_csv(data_folder+'dim_middlemen.csv')
dim_mills = pd.read_csv(data_folder+'dim_mills.csv')


dim_pickups = dim_pickups.rename(index=str, columns={"latitude": "latitude_farmer", "longitude": "longitude_farmer"})
dim_middlemen = dim_middlemen.rename(index=str, columns={"latitude": "latitude_middleman", "longitude": "longitude_middleman"})
dim_mills = dim_mills.rename(index=str, columns={"latitude": "latitude_mill", "longitude": "longitude_mill"})

# Create a dictionary of the farmers to cluster
farmer_to_cluster = dict(zip(dim_pickups['plot_id'],dim_pickups['cluster_id']))

# Calculate middleman capacity
dim_middlemen['trucks_dict']  = dim_middlemen['trucks'].map(lambda d:ast.literal_eval(d))
dim_middlemen['capacity'] = dim_middlemen['trucks_dict'].map(lambda d:np.sum([int(t)*d[t] for t in d.keys()]))


# In[7]:

# Explote data
clusters = np.unique(dim_pickups['cluster_id'])
df_clusters = pd.merge(dim_pickups,dim_pickups.pickups.apply(pd.Series).stack().reset_index(level=1, drop=True).to_frame('pickup'),left_index=True, right_index=True)
df_clusters['trucks_dict'] = df_clusters['cluster_id'].map(pd.Series(data = dim_middlemen['trucks_dict'].values, index = dim_middlemen['cluster_id']))


# In[8]:

# Calculate total capacity
dict_comparisons={}
for c in clusters:
    dict_comparisons[c] = {}
    dict_comparisons[c]['capacity'] = dim_middlemen[dim_middlemen.cluster_id == c]['capacity'].iloc[0]


# In[9]:

# Check if farmers sell more than the trucks capacity
inconsistencies = 0
quantities_changed = {}
for i, row in df_clusters.iterrows():
    max_capacity = np.max([int(k) for k in row['trucks_dict'].keys()])
    if row['volume']>max_capacity:
        df_clusters.loc[i,'volume'] = max_capacity
        quantities_changed[df_clusters.loc[i,'plot_id']] = max_capacity
        inconsistencies += 1
print("A total of %d instances where farmer produced more than truck were found and replaced" % inconsistencies,file=f)


# In[10]:

# Number of plantations picked up each day and quantities picked up each day
# Calculate the number of days
agg_quant = df_clusters.groupby(['cluster_id','pickup']).agg({'farmer_id':'count', 'volume': 'sum'})
agg_quant['overload']=agg_quant['volume']-agg_quant.apply(lambda r:dict_comparisons[r.name[0]]['capacity'],1)
outliers = (agg_quant[agg_quant['overload']>0])
print(outliers)


# In[11]:

# Create a report of inconsistent data
report = []
for index,row in outliers.iterrows():
    report.append({'cluster':index[0],
                   'farmers':row['farmer_id'],
                         'volume':row['volume'],
                         'capacity':row['volume']-row['overload'],
                        'farmer_id-plot':[r['plot_id'] for j,r in df_clusters.iterrows() if (r['cluster_id']==index[0] and r['pickup']==index[1])]})
print("Found %d trucks carrying more than their capacity" % len(outliers),file=f)
pd.DataFrame(report).to_csv(results_folder_ts+'data_cleaning_results/overcapacity.csv')
print("Report saved in data_cleaning_results/overcapacity.csv",file=f)


# In[ ]:




# In[12]:

df_clusters_c = df_clusters.copy().reset_index(drop=True)
moved_farmers = 0
# Modify inconsistent data, changing the pickup days randomly
for index,row in outliers.iterrows():
    success = False
    df_cluster_outlier = df_clusters_c[df_clusters_c['cluster_id']==index[0]]
    df_cluster_pickup_outlier = df_clusters_c[(df_clusters_c['cluster_id']==index[0])&(df_clusters_c['pickup']==index[1])]
    cluster_capacity = dict_comparisons[index[0]]['capacity']
    overload = np.sum(df_cluster_pickup_outlier['volume']) - cluster_capacity
    # Extract rows to move
    for index_to_replace, row_to_replace in df_cluster_pickup_outlier.iterrows():
        row_capacity = row_to_replace['volume']
        # Calculate the overloads of all pickup days
        df_cluster_outlier = df_clusters_c[df_clusters_c['cluster_id']==index[0]]
        agg_quant = df_cluster_outlier.groupby(['pickup']).agg({'farmer_id':'count', 'volume': 'sum'})
        agg_quant['overload']=agg_quant['volume']-cluster_capacity
        for index_to_replace_for,row_to_replace_for in agg_quant.iterrows():
            if (row_to_replace_for['overload'] + row_capacity) <= 0:
                
                df_clusters_c.loc[index_to_replace,'pickup'] = index_to_replace_for
                moved_farmers +=1
                break
        df_cluster_pickup_outlier = df_clusters_c[(df_clusters_c['cluster_id']==index[0])&(df_clusters_c['pickup']==index[1])]
        overload = np.sum(df_cluster_pickup_outlier['volume']) - cluster_capacity
        if (overload <= 0):
            success = True
            break
    if not success:
        print("Not possible to make data consistent",file=f)
print("Fixed consistency with moving schedules of %d farmers" % moved_farmers,file=f)


# In[13]:

# Parses data into numeric JSON format that can be read by C++
def json_parser(H,N,H_p,N_p,capacities,quantities,mapping, distances, n_trucks, penalties = None):
    N_ = list(range(len(N)))
    H_ = list(range(len(N),len(N)+len(H)))
    capacities_ = [int(capacities[h]*10) for h in H]
    quantities_ = [int(quantities[n]*10) for n in N]
    n_trucks_ = [int(n_trucks[h]) for h in H]
    if penalties:
        penalties_ = []
        for i,h in enumerate(H):
            penalties_.append([])
            for j,n in enumerate(N):
                penalties_[i].append(penalties[h][n])
    file_dict = {"H":H_,
                "N":N_,
                "H_p":H_p.tolist(),
                "N_p":N_p.tolist(),
                "capacities": capacities_,
                "quantities": quantities_,
                "mapping": mapping,
                "distances": distances.tolist(),
                "n_trucks": n_trucks_,
                }
    if penalties:
        file_dict.update({
                "penalties": penalties_,
            })
    return file_dict
    


# In[14]:

# Load the trucks dataset
dim_trucks = pd.read_csv(data_folder+'dim_trucks.csv')
# Change names of datasets
df_clusters_original = df_clusters.copy()
df_clusters = df_clusters_c.copy()


# In[15]:

# Number of plantations picked up each day and quantities picked up each day
# Calculate the number of days
agg_quant = df_clusters.groupby(['cluster_id','pickup']).agg({'farmer_id':'count', 'volume': 'sum'})
agg_quant['overload']=agg_quant['volume']-agg_quant.apply(lambda r:dict_comparisons[r.name[0]]['capacity'],1)
outliers = (agg_quant[agg_quant['overload']>0])
print(outliers)


# In[16]:

# Save this file into a JSON
with open(results_folder_ts+'data_cleaning_results/changed_quantities.json', 'w') as outfile:
    json.dump(quantities_changed, outfile)


# In[17]:

# Create dataset of daily production of middlemen
df_clusters.groupby(['cluster_id','pickup'])[['volume']].sum().to_csv(results_folder_ts+"middlemen_daily_production.csv")


# In[18]:

# Create a dictionary of pickups
pickups_dict = collections.defaultdict(list)
for i,row in df_clusters.iterrows():
    pickups_dict[row['plot_id']] = pickups_dict[row['plot_id']] + [row['pickup']]
with open(results_folder_ts+'data_cleaning_results/pickup_days.json', 'w') as outfile:
    json.dump(pickups_dict, outfile, indent = 4)


# In[19]:

f.close()


# In[21]:

# Check feasibility
for c in clusters:
    for d in range(days_analysis):
        # Extract mill coordinates
        M_p = np.array(dim_mills[dim_mills['code'] == 'SKIP'][["latitude_mill","longitude_mill"]].values)
        trucks = dim_trucks[dim_trucks['cluster_id']==c]['truck_id'].values
        k = len(trucks)
        H = ['h_'+str(i) for i in range(k)]
        H_p = np.array(dim_middlemen[dim_middlemen['cluster_id'] == c][["latitude_middleman","longitude_middleman"]])
        H_p = np.tile(H_p,(k,1))
        capacities = dict(zip(H,dim_trucks[dim_trucks['cluster_id']==c]['capacity'].values))
        mapping_trucks = dict(zip(H,trucks))
        df_cluster_day = df_clusters[(df_clusters['pickup'] == d) & (df_clusters['cluster_id'] == c)]
        n = len(df_cluster_day)
        m = 1
        if (n!=0):
            N = ['n_'+str(i) for i in range(n)]
            M = ['m_'+str(i) for i in range(m)]
            quantities = {f: df_cluster_day["volume"].values[i] for i,f in enumerate(N)}
            N_p = np.array(df_cluster_day[['latitude_farmer','longitude_farmer']])
            # Create VRP
            vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, 'road', M, M_p, distance_matrix = np.zeros([len(N)+len(H),len(N)+len(H)]))
            assert( vrp.feasibility())


# In[22]:

# Create daily instances
bash_path = results_folder_ts+"bash/bash_daily.sh"
open(bash_path, 'w').close()
for c in clusters:
    for d in range(days_analysis):
        # Extract mill coordinates
        M_p = np.array(dim_mills[dim_mills['code'] == 'SKIP'][["latitude_mill","longitude_mill"]].values)
        trucks = dim_trucks[dim_trucks['cluster_id']==c]['truck_id'].values
        k = len(trucks)
        H = ['h_'+str(i) for i in range(k)]
        H_p = np.array(dim_middlemen[dim_middlemen['cluster_id'] == c][["latitude_middleman","longitude_middleman"]])
        H_p = np.tile(H_p,(k,1))
        capacities = dict(zip(H,dim_trucks[dim_trucks['cluster_id']==c]['capacity'].values))
        mapping_trucks = dict(zip(H,trucks))
        df_cluster_day = df_clusters[(df_clusters['pickup'] == d) & (df_clusters['cluster_id'] == c)]
        n = len(df_cluster_day)
        m = 1
        if (n!=0):
            N = ['n_'+str(i) for i in range(n)]
            M = ['m_'+str(i) for i in range(m)]
            quantities = {f: df_cluster_day["volume"].values[i] for i,f in enumerate(N)}
            N_p = np.array(df_cluster_day[['latitude_farmer','longitude_farmer']])
            # Create VRP
            vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, 'road', M, M_p)
            
            # Dump it in case we want to reuse it
            pickle.dump(vrp, open( results_folder_ts+"vrps/daily/daily_"+"cluster_"+str(c)+"_day_"+str(d)+".p", "wb" ) )

            mapping = {h:(mapping_trucks[h],d, c) for h in H}
            for i,ind in enumerate(df_cluster_day['plot_id']):
                mapping[N[i]] = ind
            n_trucks = dict(zip(H,[1]*len(H)))
            file_dict = json_parser(H,N,H_p,N_p,capacities,quantities,mapping, vrp.distance, n_trucks)
            with open( results_folder_ts+"instances/daily/daily_"+"cluster_"+str(c)+"_day_"+str(d)+".json", "wb" ) as fp:
                json.dump(file_dict, fp)
            file_dict_pickle = {
                "H":H,
                "N":N,
                "M":M,
                "H_p":H_p,
                "N_p":N_p,
                "M_p":M_p,
                "capacities":capacities,
                "quantities":quantities,
                "mapping":mapping,
                "distance":vrp.distance,
                "type_dist":"road",
                "n_trucks":n_trucks,
            }
            
            if not vrp.feasibility():
                print("Not feasible")
                print(quantities)
                print(capacities)
                print (mapping)
            
            pickle.dump( file_dict_pickle, open( results_folder_ts+"instances/daily/daily_"+"cluster_"+str(c)+"_day_"+str(d)+".p", "wb" ) )
            with open(bash_path, 'a') as file:
                file.write("python2 "+folder_project+"Code/PipelineScripts/load_and_solve_instance.py "+"instances/daily/daily_"+"cluster_"+str(c)+"_day_"+str(d) + " "+ folder_project +" Results/"+ts+"/"+"\n")
                


# In[23]:

# Create temporal instances
bash_path = results_folder_ts+"bash/bash_temporal.sh"
open(bash_path, 'w').close()

for c in clusters:
    #print("Processing config files of cluster: %d" % (c),file=f)
    M_p = np.array(dim_mills[dim_mills['code'] == 'SKIP'][["latitude_mill","longitude_mill"]].values)
    trucks = dim_trucks[dim_trucks['cluster_id']==c]['truck_id'].values
    trucks
    k = len(trucks)
    H = ['h_'+str(i) for i in range(k)]
    H_p = np.array(dim_middlemen[dim_middlemen['cluster_id'] == c][["latitude_middleman","longitude_middleman"]])
    H_p = np.tile(H_p,(k,1))
    H_p
    capacities = dict(zip(H,dim_trucks[dim_trucks['cluster_id']==c]['capacity'].values))
    capacities
    mapping_trucks = dict(zip(H,trucks))
    df_cluster = df_clusters[(df_clusters['cluster_id'] == c)]

    # Construct optimization problem
    n = len(df_cluster)
    m = 1

    N = ['n_'+str(i) for i in range(n)]
    M = ['m_'+str(i) for i in range(m)]


    # Construct the mapping
    mapping_trucks = dict(zip(H,trucks))
    mapping = {h:(mapping_trucks[h],c) for h in H}
    for i,ind in enumerate(df_cluster['plot_id']):
        mapping[N[i]] = ind

    quantities = {f: df_cluster['volume'].values[i] for i,f in enumerate(N)}
    N_p = np.array(df_cluster[['latitude_farmer','longitude_farmer']])
    vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, 'road', M, M_p)
    # Dump it in case we want to reuse it
    pickle.dump(vrp, open( results_folder_ts+"vrps/temporal/temporal_"+"cluster_"+str(c)+".p", "wb" ) )

    n_trucks = dict(zip(H,[days_analysis]*len(H)))
    file_dict = json_parser(H,N,H_p,N_p,capacities,quantities,mapping, vrp.distance, n_trucks)
    with open( results_folder_ts+"instances/temporal/temporal_"+"cluster_"+str(c)+".json", "wb" ) as fp:
        json.dump(file_dict, fp)
    file_dict_pickle = {
        "H":H,
        "N":N,
        "M":M,
        "H_p":H_p,
        "N_p":N_p,
        "M_p":M_p,
        "capacities":capacities,
        "quantities":quantities,
        "mapping":mapping,
        "distance":vrp.distance,
        "type_dist":"road",
        "n_trucks":n_trucks
    }
    pickle.dump( file_dict_pickle, open( results_folder_ts+"instances/temporal/temporal_"+"cluster_"+str(c)+".p", "wb" ) )
    with open(bash_path, 'a') as file:
        file.write("python2 "+folder_project+"Code/PipelineScripts/load_and_solve_instance.py "+"instances/temporal/temporal_"+"cluster_"+str(c) + " "+ folder_project +" Results/"+ts+"/"+"\n")


# In[24]:

# We now solve the spatial problem
# Results spatial
bash_path = results_folder_ts+"bash/bash_spatial.sh"
open(bash_path, 'w').close()

# Extract all trucks that are available
capacities = {}
H_p = np.empty([0,2])
H = []; i = 0
truck_mapping = {}
for c in clusters:
    middle_df_c = dim_middlemen[dim_middlemen.cluster_id == c]
    trucks = dim_trucks[dim_trucks['cluster_id']==c]['truck_id'].values
    capacities_v = dim_trucks[dim_trucks['cluster_id']==c]['capacity'].values
    for j,t in enumerate(trucks):
        H_p = np.vstack([H_p, np.array(middle_df_c[['latitude_middleman','longitude_middleman']])])
        H.append('h_'+str(i))
        truck_mapping['h_'+str(i)] = (t,c)
        capacities['h_'+str(i)] = capacities_v[j]
        i += 1
k = len(capacities)
M_p = np.array(dim_mills[dim_mills['code'] == 'SKIP'][["latitude_mill","longitude_mill"]].values)
for d in range(days_analysis):
    #print("Processing config files of dat: %d" % (d),file=f)
    df_day = df_clusters[(df_clusters['pickup'] == d)]
    print(len(df_day))
    n = len(df_day)
    m = 1
    N = ['n_'+str(i) for i in range(n)]
    M = ['m_'+str(i) for i in range(m)]
    mapping = {}
    for h in H:
        mapping[h] = (truck_mapping[h][0],d,truck_mapping[h][1])
    for i,ind in enumerate(df_day['plot_id']):
        mapping[N[i]] = ind
        
    # Construct optimization problem
    quantities = {f: df_day['volume'].values[i] for i,f in enumerate(N)}
    N_p = np.array(df_day[['latitude_farmer','longitude_farmer']])
    vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, 'road', M, M_p)
    
    # Create the penalty vector for matching middlemen with farmers that do not correspond to them
    penalties = {}
    for h in H:
        penalties[h] = {}
        for n in N:
            cluster_n = farmer_to_cluster[mapping[n]]
            if truck_mapping[h][1] == cluster_n:
                penalties[h][n] = 0
            else:
                penalties[h][n] = quantities[n]
    
    
    # Dump it in case we want to reuse it
    pickle.dump(vrp, open( results_folder_ts+"vrps/spatial/spatial_"+"day_"+str(d)+".p", "wb" ) )

    n_trucks = dict(zip(H,[1]*len(H)))

    file_dict = json_parser(H,N,H_p,N_p,capacities,quantities,mapping, vrp.distance,n_trucks, penalties)
    #print(file_dict)

    with open( results_folder_ts+"instances/spatial/spatial_"+"day_"+str(d)+".json", "wb" ) as fp:
        json.dump(file_dict, fp)

    file_dict_pickle = {
        "H":H,
        "N":N,
        "M":M,
        "H_p":H_p,
        "N_p":N_p,
        "M_p":M_p,
        "capacities":capacities,
        "quantities":quantities,
        "mapping":mapping,
        "distance":vrp.distance,
        "type_dist":"road",
        "n_trucks":n_trucks,
        "penalties":penalties,
    }

    pickle.dump( file_dict_pickle, open( results_folder_ts+"instances/spatial/spatial_"+"day_"+str(d)+".p", "wb" ) )

    with open(bash_path, 'a') as file:
        file.write("python2 "+folder_project+"Code/PipelineScripts/load_and_solve_instance.py "+"instances/spatial/spatial_"+"day_"+str(d) + " "+ folder_project +" Results/"+ts+"/"+"\n")
        


# In[29]:

# Create bash scripts for the cluster

penalty = float(sys.argv[2])

header = """#!/bin/bash
#
#SBATCH --job-name=test
#
#SBATCH --time=%s:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G

cd ~/vrp_project
source ~/.virtualenvs/vrp_project/bin/activate
module load gurobi/8.0.1_py27
srun python ~/vrp_project/Code/PipelineScripts/load_and_solve_instance.py instances/%s/%s ~/vrp_project/ %s %d -p normal,dev,gpu
"""


folder_path = results_folder_ts + 'bash/'
cluster_path = "~/vrp_project/" + results_folder_ts + 'bash/cluster/'


submission = folder_path + 'submission_daily.sh'
open(submission, 'w').close()

for c in clusters:
    for d in range(days_analysis):
        instance = "daily_"+"cluster_"+str(c)+"_day_"+str(d)
        path = folder_path + '/cluster/' + instance + '.sh'
        open(path, 'w').close()
        with open(path, 'a') as file:
            file.write(header % ('00:20',"daily",instance,results_folder_ts,penalty))
        with open(submission, 'a') as file:
            file.write('sbatch ' + cluster_path + instance + '.sh' + '\n')

submission = folder_path + 'submission_temporal.sh'
open(submission, 'w').close()
for c in clusters:
    instance = "temporal_"+"cluster_"+str(c)
    path = folder_path + '/cluster/' + instance + '.sh'
    open(path, 'w').close()
    with open(path, 'a') as file:
        file.write(header % ('01:30',"temporal",instance,results_folder_ts,penalty))
    with open(submission, 'a') as file:
        file.write('sbatch ' + cluster_path + instance + '.sh' + '\n')
        
submission = folder_path + 'submission_spatial.sh'
open(submission, 'w').close()
for d in range(days_analysis):
    instance = "spatial_"+"day_"+str(d)
    path = folder_path + '/cluster/' + instance + '.sh'
    open(path, 'w').close()
    with open(path, 'a') as file:
        file.write(header % ('01:30',"spatial",instance,results_folder_ts,penalty))
    with open(submission, 'a') as file:
        file.write('sbatch ' + cluster_path + instance + '.sh' + '\n')
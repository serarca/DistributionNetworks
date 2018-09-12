from __future__ import print_function
import sys

# Get season from argument
season = sys.argv[1]
print("Season is: %s" % season)
with open(('report_%s.txt' % season), 'w') as f:

    # coding: utf-8

    # In[33]:


    #%matplotlib inline
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


    sys.path.insert(0, '/Users/sergiocamelo/Dropbox/Sergio-Joann/Code')
    sys.path.insert(0, '/Users/sergiocamelo/Dropbox/Sergio-Joann/Code/LowerBoundsC')

    import solver as solver
    import distances as distances
    import VRPClass


    # In[34]:


    data_folder = '/Users/sergiocamelo/Dropbox/Sergio-Joann/StandardizedData/'
    period = 3*28 # Use three month period
    days_analysis = 14


    # In[35]:


    season_volume = season + '_volume'


    # In[36]:


    # Load datasets
    dim_farmers = pd.read_csv(data_folder+'dim_farmers.csv')
    dim_middlemen = pd.read_csv(data_folder+'dim_middlemen.csv')
    dim_mills = pd.read_csv(data_folder+'dim_mills.csv')
    harvest_frequency_mapping = pd.read_csv(data_folder+'harvest_frequency_mapping.csv')
    print("Number of plantations: %d" % (len(dim_farmers)),file=f)
    print("Number of middlemen: %d" % (len(dim_middlemen)),file=f)


    # In[5]:


    # Create a unique identifier for the farmer plot
    dim_farmers['plot_id'] = [dim_farmers['farmer_id'][i] + '-'+str(dim_farmers['plot_number'][i]) for i in range(len(dim_farmers))]


    # In[6]:


    # Calculate middleman capacity
    dim_middlemen['trucks_dict']  = dim_middlemen['trucks'].map(lambda d:ast.literal_eval(d))
    dim_middlemen['capacity'] = dim_middlemen['trucks_dict'].map(lambda d:np.sum([int(t)*d[t] for t in d.keys()]))


    # In[7]:


    # Join farmers and middlemen data
    # Change names of columns
    dim_farmers = dim_farmers.rename(index=str, columns={"latitude": "latitude_farmer", "longitude": "longitude_farmer"})
    dim_middlemen = dim_middlemen.rename(index=str, columns={"latitude": "latitude_middleman", "longitude": "longitude_middleman"})
    dim_mills = dim_mills.rename(index=str, columns={"latitude": "latitude_mill", "longitude": "longitude_mill"})

    result = pd.merge(dim_farmers, dim_middlemen, on=['cluster_id'], how='inner')
    print("Number of plantations with middleman: %d" % (len(result)),file=f)


    # In[8]:


    # Check if there are any duplicates
    result[result.duplicated(subset=['farmer_id','plot_number'], keep=False)].to_csv(data_folder+'data_cleaning_results/duplicates.csv')
    print("Total of duplicates: %d" % (len(result[result.duplicated(subset=['farmer_id','plot_number'], keep=False)])),file=f)
    print ("Duplicates saved in data_cleaning_results/duplicates.csv",file=f)


    # In[9]:


    #Use only data with lat_lon and with productions
    df_full = result[np.logical_and(pd.notnull(result['longitude_farmer']),pd.notnull(result['latitude_farmer']))].copy()
    print("Number of plantations with latlon: %d" % (len(df_full)),file=f)
    # Use data with productions
    df_full = df_full[df_full[season+'_rate']!=0].copy()
    df_full = df_full[pd.notnull(df_full[season+'_rate'])].copy()
    df_full = df_full[df_full[season+'_volume']!=0].copy()
    print("Number of plantations that produce during the season : %d" % (len(df_full)),file=f)


    # In[10]:


    # Map number of days
    harvest_frequency_mapping_dict = {row[0]:row[1] for i,row in harvest_frequency_mapping.iterrows()}
    df_full['rate'] = df_full[season+'_rate'].map(harvest_frequency_mapping_dict)


    # In[11]:


    # Has pickup date 
    df_full = df_full[pd.notnull(df_full['date_last_sold'])].copy()
    print("Number of plantations with last date: %d" % (len(df_full)),file=f)


    # In[12]:


    # Generate days of pickup
    ref_day = datetime.datetime.strptime('1/1/2000', "%m/%d/%Y")
    days = np.array([(parse(v)-ref_day).days for v in df_full['date_last_sold'].values])
    df_full['day_mod'] = days%period
    def calculate_pickup_days(row):
        d = row['day_mod']
        freq = row[season+'_rate']
        l = []
        for i in range(int(period/freq)):
            l.append((d + i * freq)%period)
        return l
    df_full['pickup_days'] = df_full.apply(calculate_pickup_days, axis=1)


    # In[13]:


    # Explote data
    clusters = np.unique(df_full['cluster_id'])
    df_exploted = pd.merge(df_full,df_full.pickup_days.apply(pd.Series).stack().reset_index(level=1, drop=True).to_frame('pickup'),left_index=True, right_index=True)
    for c in clusters:
        print("Number of plantations of cluster %d: %d" % (c,len(df_full[df_full['cluster_id'] == c].copy())),file=f)
    df_clusters = df_exploted[(np.array([c in clusters for c in df_exploted['cluster_id']])) & (df_exploted['pickup'] < days_analysis)].copy()


    # In[14]:


    # Print middlemen data in the results folder
    dim_middlemen[[dim_middlemen['cluster_id'][i] in clusters for i in range(len(dim_middlemen))]][['cluster_id','trucks','mills','capacity']].to_csv(data_folder+"tables_for_report/middlemen_data.csv", index=False)


    # In[15]:


    # Calculate total capacity
    dict_comparisons={}
    for c in clusters:
        dict_comparisons[c] = {}
        dict_comparisons[c]['capacity'] = dim_middlemen[dim_middlemen.cluster_id == c]['capacity'].iloc[0]


    # In[16]:


    # Round producing quantities to the decimal up
    df_clusters[season + '_volume'] = np.ceil(df_clusters[season + '_volume'] * 10)/10


    # In[17]:


    # Number of plantations picked up each day and quantities picked up each day
    # Calculate the number of days
    agg_quant = df_clusters.groupby(['cluster_id','pickup']).agg({'farmer_id':'count', season+'_volume': 'sum'})
    agg_quant['overload']=agg_quant[season+'_volume']-agg_quant.apply(lambda r:dict_comparisons[r.name[0]]['capacity'],1)
    outliers = (agg_quant[agg_quant['overload']>0])
    print(outliers)


    # In[18]:


    # Create a report of inconsistent data
    report = []
    for index,row in outliers.iterrows():
        report.append({'cluster':index[0],
                       'farmers':row['farmer_id'],
                             season+'_volume':row[season+'_volume'],
                             'capacity':row[season+'_volume']-row['overload'],
                            'farmer_id-plot':[r['farmer_id']+'-'+str(r['plot_number']) for j,r in df_clusters.iterrows() if (r['cluster_id']==index[0] and r['pickup']==index[1])]})
    print("Found %d trucks carrying more than their capacity" % len(outliers),file=f)
    pd.DataFrame(report).to_csv(data_folder+'data_cleaning_results/overcapacity.csv')
    print("Report saved in data_cleaning_results/overcapacity.csv",file=f)


    # In[19]:


    df_clusters_c = df_clusters.copy()
    moved_farmers = 0
    # Modify inconsistent data, changing the pickup days randomly
    for index,row in outliers.iterrows():
        success = False
        df_cluster_outlier = df_clusters_c[df_clusters_c['cluster_id']==index[0]]
        df_cluster_pickup_outlier = df_clusters_c[(df_clusters_c['cluster_id']==index[0])&(df_clusters_c['pickup']==index[1])]
        cluster_capacity = dict_comparisons[index[0]]['capacity']
        overload = np.sum(df_cluster_pickup_outlier[season+'_volume']) - cluster_capacity
        # Extract rows to move
        for index_to_replace, row_to_replace in df_cluster_pickup_outlier.iterrows():
            row_capacity = row_to_replace[season+'_volume']
            # Calculate the overloads of all pickup days
            df_cluster_outlier = df_clusters_c[df_clusters_c['cluster_id']==index[0]]
            agg_quant = df_cluster_outlier.groupby(['pickup']).agg({'farmer_id':'count', season+'_volume': 'sum'})
            agg_quant['overload']=agg_quant[season+'_volume']-cluster_capacity
            for index_to_replace_for,row_to_replace_for in agg_quant.iterrows():
                if (row_to_replace_for['overload'] + row_capacity) <= 0:
                    df_clusters_c.loc[index_to_replace,'pickup'] = index_to_replace_for
                    moved_farmers +=1
                    break
            df_cluster_pickup_outlier = df_clusters_c[(df_clusters_c['cluster_id']==index[0])&(df_clusters_c['pickup']==index[1])]
            overload = np.sum(df_cluster_pickup_outlier[season+'_volume']) - cluster_capacity
            if (overload <= 0):
                success = True
                break
        if not success:
            print("Not possible to make data consistent",file=f)
    print("Fixed consistency with moving schedules of %d farmers" % moved_farmers,file=f)


    # In[20]:


    # Parses data into numeric JSON format that can be read by C++
    def json_parser(H,N,H_p,N_p,capacities,quantities,mapping, distances, n_trucks):
        N_ = list(range(len(N)))
        H_ = list(range(len(N),len(N)+len(H)))
        capacities_ = [int(capacities[h]*10) for h in H]
        quantities_ = [int(quantities[n]*10) for n in N]
        n_trucks_ = [int(n_trucks[h]) for h in H]
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
        return file_dict
        


    # In[21]:


    # # Create trucks dataset (Only run once)
    # truck_dicts = []
    # j = 0
    # for i,r in dim_middlemen.iterrows():
    #     truck_dic = ast.literal_eval(r['trucks'])
    #     for capacity in truck_dic.keys():
    #         for i in range(truck_dic[capacity]):
    #             truck_dicts.append({
    #                     "cluster_id":r['cluster_id'],
    #                     "truck_id":"t_"+str(j),
    #                     "capacity":int(capacity)
    #                 })
    #             j += 1
    # pd.DataFrame(truck_dicts).to_csv(data_folder+"dim_trucks.csv",index=False)


    # In[22]:


    # Load the trucks dataset
    dim_trucks = pd.read_csv(data_folder+'dim_trucks.csv')
    # Change names of datasets
    df_clusters_oiginal = df_clusters.copy()
    df_clusters = df_clusters_c.copy()


    # In[23]:


    # Check if farmers sell more than the trucks capacity
    inconsistencies = 0
    for i, row in df_clusters.iterrows():
        max_capacity = np.max([int(k) for k in row['trucks_dict'].keys()])
        if row[season+'_volume']>max_capacity:
            df_clusters.loc[i,season+'_volume'] = max_capacity
            inconsistencies += 1
    print("A total of %d instances where farmer produced more than truck were found and replaced" % inconsistencies,file=f)


    # In[24]:


    # Create daily instances
    bash_path = data_folder+"instances/bash_daily.sh"
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
                quantities = {f: df_cluster_day[season+"_volume"].values[i] for i,f in enumerate(N)}
                N_p = np.array(df_cluster_day[['latitude_farmer','longitude_farmer']])
                vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, 'road', M, M_p)
                mapping = {h:(mapping_trucks[h],d, c) for h in H}
                for i,ind in enumerate(df_cluster_day['plot_id']):
                    mapping[N[i]] = ind
                n_trucks = dict(zip(H,[1]*len(H)))
                file_dict = json_parser(H,N,H_p,N_p,capacities,quantities,mapping, vrp.distance, n_trucks)
                with open( data_folder+"instances/daily/daily_"+"cluster_"+str(c)+"_day_"+str(d)+".json", "wb" ) as fp:
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
                pickle.dump( file_dict_pickle, open( data_folder+"instances/daily/daily_"+"cluster_"+str(c)+"_day_"+str(d)+".p", "wb" ) )
                with open(bash_path, 'a') as file:
                    file.write("python2 /Users/sergiocamelo/Dropbox/Sergio-Joann/Code/DataInstances/load_and_solve_instance.py "+data_folder+"instances/daily/daily_"+"cluster_"+str(c)+"_day_"+str(d) + "\n")


    # In[25]:


    # Create temporal instances
    bash_path = data_folder+"instances/bash_temporal.sh"
    open(bash_path, 'w').close()

    for c in clusters:
        print("Processing config files of cluster: %d" % (c),file=f)
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

        quantities = {f: df_cluster[season_volume].values[i] for i,f in enumerate(N)}
        N_p = np.array(df_cluster[['latitude_farmer','longitude_farmer']])
        vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, 'road', M, M_p)
        n_trucks = dict(zip(H,[days_analysis]*len(H)))
        file_dict = json_parser(H,N,H_p,N_p,capacities,quantities,mapping, vrp.distance, n_trucks)
        with open( data_folder+"instances/temporal/temporal_"+"cluster_"+str(c)+".json", "wb" ) as fp:
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
        pickle.dump( file_dict_pickle, open( data_folder+"instances/temporal/temporal_"+"cluster_"+str(c)+".p", "wb" ) )
        with open(bash_path, 'a') as file:
            file.write("python2 /Users/sergiocamelo/Dropbox/Sergio-Joann/Code/DataInstances/load_and_solve_instance.py "+data_folder+"instances/temporal/temporal_"+"cluster_"+str(c) + "\n")


    # In[26]:


    # We now solve the spatial problem
    # Results spatial
    bash_path = data_folder+"instances/bash_spatial.sh"
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
        print("Processing config files of dat: %d" % (d),file=f)
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
        quantities = {f: df_day[season_volume].values[i] for i,f in enumerate(N)}
        N_p = np.array(df_day[['latitude_farmer','longitude_farmer']])
        vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, 'road', M, M_p)
        n_trucks = dict(zip(H,[1]*len(H)))

        file_dict = json_parser(H,N,H_p,N_p,capacities,quantities,mapping, vrp.distance,n_trucks)


        with open( data_folder+"instances/spatial/spatial_"+"day_"+str(d)+".json", "wb" ) as fp:
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

        pickle.dump( file_dict_pickle, open( data_folder+"instances/spatial/spatial_"+"day_"+str(d)+".p", "wb" ) )

        with open(bash_path, 'a') as file:
            file.write("python2 /Users/sergiocamelo/Dropbox/Sergio-Joann/Code/DataInstances/load_and_solve_instance.py "+data_folder+"instances/spatial/spatial_"+"day_"+str(d) + "\n")


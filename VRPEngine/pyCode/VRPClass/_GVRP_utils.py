import numpy as np
import road_distance
import networkx as nx

def coordinates(nodes):
    locations = []
    node_ids = []
    for n in nodes:
        locations.append(list(n.location))
        node_ids.append(n.node_id)
    return {
        "locs":np.array(locations),
        "node_ids":node_ids
    }

def f_penalty(truck, farmer, penalty_type):
    if penalty_type == "village":
        if (farmer.village == truck.village) or (farmer.cluster_id_sq == truck.cluster_id):
            return 0.0
        else:
            return 1.0

# Price of the load of a truck at a mill asssuming truck is 90% loaded
def load_price(mill, truck):
    return truck.capacity * 1000.0 * mill.price * 0.9


@classmethod
def create_VRP_instance(cls, ms,fs,ts,metadata,dist_type="road",network = None, penalty_type = None, truck_types = ["2_ton_diesel","9_ton"], weather_conditions="norain"):
    # Process mills
    M_coord = coordinates(ms)
    M_p = M_coord['locs']
    M_ids =  M_coord['node_ids']
    M = ['m_'+str(i) for i in range(len(M_p))]
    mills_dict = {}
    for m in ms:
        mills_dict[m.node_id] = m

    # Process farmers
    N_coord = coordinates(fs)
    N_p = N_coord['locs']
    N_ids =  N_coord['node_ids']
    N = ['n_'+str(i) for i in range(len(N_p))]
    # Get their quantities
    qs = {'n_'+str(i):f.quantity for i,f in enumerate(fs)}


    # Process trucks
    H_coord = coordinates(ts)
    H_p = H_coord['locs']
    H_ids =  H_coord['node_ids']
    H = ['h_'+str(i) for i in range(len(H_p))]
    # Get their capacities
    cs = {'h_'+str(i):t.capacity for i,t in enumerate(ts)}
    n_trucks = dict(zip(H,[1]*len(H)))

    # Construct closest mill
    closest_mill = {}
    for f in fs:
        kf = f.node_id
        closest_mill[kf] = {}
        for t in ts:
            kt = t.node_id
            min_dist = np.float('inf')
            min_mill = ''
            for m in ms:
                km = m.node_id
                if dist_type=="road_dist":
                    dist = road_distance.road_distance(f.location,m.location)['distance'] + road_distance.road_distance(m.location,t.location)['distance']
                elif dist_type=="road_time":
                    dist = road_distance.road_distance(f.location,m.location)['time'] + road_distance.road_distance(m.location,t.location)['time']
                elif dist_type=="network_dist":
                    weight = "weight_dist"
                    dist = network.network_distance(f.location,m.location,weight) + network.network_distance(m.location,t.location,weight)
                elif dist_type=="network_time":
                    weight = "weight_time"
                    dist = network.network_distance(f.location,m.location,weight) + network.network_distance(m.location,t.location,weight)
                elif dist_type=="zero_dist":
                    dist = 0
                ## This should be changed at some point
                elif dist_type=="network_cost":
                    weight = t.truck_type + "_" + weather_conditions
                    dist = network.network_distance(f.location,m.location,weight)  + network.network_distance(m.location,t.location,weight)
                else:
                    assert(False)
                if dist<min_dist:
                    min_dist = dist
                    min_mill = km
            closest_mill[kf][kt] = min_mill

    # Create mapping
    mapping = {
        "mills":M_ids,
        "trucks":H_ids,
        "farmers":N_ids,
        "closest_mill":closest_mill,
        "metadata":metadata
    }

    # Create penalty
    penalties = None
    if penalty_type:
        penalties = {}
        for i,truck in enumerate(ts):
            penalties['h_'+str(i)] = {}
            for j,farmer in enumerate(fs):
                penalties['h_'+str(i)]['n_'+str(j)] = f_penalty(truck, farmer, penalty_type)


    # Create distances
    nodes_list = fs + ts
    network_distance_matrix = np.zeros([len(nodes_list),len(nodes_list)])
    if dist_type == 'network_cost':
        network_distance_matrices = {}
        for truck_type in truck_types:
            weight = truck_type + "_" +weather_conditions
            network_distance_matrices[truck_type] = np.zeros([len(nodes_list),len(nodes_list)])
            for i1,n1 in enumerate(nodes_list):
                for i2,n2 in enumerate(nodes_list):
                    if i1!=i2:
                        if (n1.t == "farmer" and n2.t == "truck"):
                            #import pdb; pdb.set_trace()
                            network_distance_matrices[truck_type][i1,i2] = (
                                network.network_distance(n1.location,mills_dict[closest_mill[n1.node_id][n2.node_id]].location,weight)+
                                network.network_distance(mills_dict[closest_mill[n1.node_id][n2.node_id]].location,n2.location,weight)+
                                n2.fixed_cost/2.0
                            )
                        elif (n1.t == "truck" and n2.t == "farmer"):
                             network_distance_matrices[truck_type][i1,i2] = (
                                network.network_distance(n1.location,n2.location,truck_type)+
                                n1.fixed_cost/2.0
                            )
                        else:
                            network_distance_matrices[truck_type][i1,i2] = (
                                network.network_distance(n1.location,n2.location,truck_type)
                            )
        vrp = cls(H, N, H_p, N_p, qs, cs, dist_type, 
                   M=M, M_p=M_p, distance_matrix = network_distance_matrices, n_trucks = n_trucks, mapping = mapping, network=network, penalties = penalties)




    elif dist_type == 'network_dist' or dist_type == 'network_time':
        if dist_type == 'network_dist':
            weight = "weight_dist"
        elif dist_type == 'network_time':
            weight = "weight_time"
        for i1,n1 in enumerate(nodes_list):
            for i2,n2 in enumerate(nodes_list):
                if i1!=i2:
                    if (n1.t == "farmer" and n2.t == "truck"):
                        #import pdb; pdb.set_trace()
                        network_distance_matrix[i1,i2] = (
                            network.network_distance(n1.location,mills_dict[closest_mill[n1.node_id][n2.node_id]].location,weight)+
                            network.network_distance(mills_dict[closest_mill[n1.node_id][n2.node_id]].location,n2.location,weight)
                        )
                    else:
                        network_distance_matrix[i1,i2] = (
                            network.network_distance(n1.location,n2.location,weight)
                        )
        vrp = cls(H, N, H_p, N_p, qs, cs, dist_type, 
                   M=M, M_p=M_p, distance_matrix = network_distance_matrix, n_trucks = n_trucks, mapping = mapping, network=network, penalties = penalties)

    elif dist_type == "road_dist":
        vrp = cls(H, N, H_p, N_p, qs, cs, 'road_dist', 
                   M=M, M_p=M_p, n_trucks = n_trucks, mapping = mapping, penalties = penalties)
    elif dist_type == "road_time":
        vrp = cls(H, N, H_p, N_p, qs, cs, 'road_time', 
                   M=M, M_p=M_p, n_trucks = n_trucks, mapping = mapping, penalties = penalties)
    elif dist_type == "zero_dist":
        vrp = cls(H, N, H_p, N_p, qs, cs, 'zero_dist', 
                   M=M, M_p=M_p, n_trucks = n_trucks, mapping = mapping, penalties = penalties)
    
    return vrp
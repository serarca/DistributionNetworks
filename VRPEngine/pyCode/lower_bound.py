# In this code we implement the lower bounds of Baldacci's implementation
# Import VRP package
import sys
sys.path.insert(0, '/Users/sergiocamelo/Dropbox/Sergio-Joann/Code')
import VRPClass
import ipdb

import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import copy
from copy import deepcopy
import cpp_lower_bounds
import ipdb
from gurobipy import *
import copy



def possible_values(quantities, maximum):
    possible_values = []
    q_int = [int(quantities[n]) for n in quantities.keys()]
    def GCD(a, b):
        if b == 0:
            return a
        else:
            return GCD(b, a % b)
    reducer = int(reduce(GCD, (q_int)))
    s = reducer
    while s<=maximum:
        possible_values.append(s)
        s+=reducer
    return possible_values

# Generates a dictionary of reduced costs given a distance, a lambda vector
def reduced_cost_dict(lamb, distance_dictionary, N):
    distance = copy.deepcopy(distance_dictionary)
    for k1 in distance.keys():
        for k2 in distance[k1].keys():
            if k1 in N:
                distance[k1][k2] = distance[k1][k2] - 1.0/2*lamb[k1]
            if k2 in N:
                distance[k1][k2] = distance[k1][k2] - 1.0/2*lamb[k2]
    return distance

# Generates a dictionary of reduced costs given a distance, a lambda vector
def reduced_cost_dict_complete(lamb, mu, distance_dictionary, N, H):
    distance = copy.deepcopy(distance_dictionary)
    for k1 in distance.keys():
        for k2 in distance[k1].keys():
            if k1 in N:
                distance[k1][k2] = distance[k1][k2] - 1.0/2*lamb[k1]
            elif k1 in H:
                distance[k1][k2] = distance[k1][k2] - 1.0/2*mu[k1]
            if k2 in N:
                distance[k1][k2] = distance[k1][k2] - 1.0/2*lamb[k2]
            elif k2 in H:
                distance[k1][k2] = distance[k1][k2] - 1.0/2*mu[k2]
    return distance

# Generates a matrix of reduced costs given a distance and a lambda vector
def reduced_cost_mat(lamb_, distance_mat, N_):
    distance = copy.deepcopy(distance_mat)
    for i in range(len(N_)):
        for j in range(distance_mat.shape[0]):
            distance[i][j] -= 1.0/2*lamb_[i]
    for i in range(len(N_)):
        for j in range(distance_mat.shape[0]):
            distance[j][i] -= 1.0/2*lamb_[i]
    return distance.astype('float64')

# Generates a matrix of reduced costs given a distance and a lambda vector and a mu vector
def reduced_cost_mat_complete(distance_mat, lamb_, mu_):
    penalties = np.concatenate((lamb_, mu_))
    distance = copy.deepcopy(distance_mat)
    for i in range(len(penalties)):
        for j in range(len(penalties)):
            distance[i][j] -= (1.0/2*penalties[i] + 1.0/2*penalties[j])
    return distance.astype('float64')

def construct_q_paths(h,truck_capacity,N,distance,values,values_pos,quantities,direction):

    # Initialize the routes
    f = {}; phi = {}; p = {}; q_route = {}; q_route_2 = {};
    for l in range(len(values)):
        f[l] = {}; phi[l] = {}; p[l] = {}; q_route[l] = {}; q_route_2[l] = {};
        for n in N:
            f[l][n] = float('inf'); phi[l][n] = float('inf'); p[l][n] = float('inf');
            q_route[l][n] = []; q_route_2[l][n] = [];

    # Initialize the routes
    for n in N:
        q = quantities[n]
        if q<= truck_capacity:
            l = values_pos[q]
            if direction == 'left':
                f[l][n] = distance[h][n]
            else:
                f[l][n] = distance[n][h]
            p[l][n] = h
            q_route[l][n] = [h, n]

    # Calculate the recursion
    for l in range(len(values)):
        Q = values[l]
        g = {}
        g_type = {}
        for x_i in N:
            q_p = Q - quantities[x_i]
            g[x_i] = {}
            g_type[x_i] = {}
            if q_p>0:
                l_p = values_pos[q_p]
                for x_j in N:
                    if x_i!=x_j:
                        if p[l_p][x_j]!=x_i:
                            g[x_i][x_j] = f[l_p][x_j] + distance[x_j][x_i]
                            g_type[x_i][x_j] = q_route[l_p][x_j] + [x_i]
                        else:
                            g[x_i][x_j] = phi[l_p][x_j] + distance[x_j][x_i]
                            g_type[x_i][x_j] = q_route_2[l_p][x_j] + [x_i]
        for x_i in N:
            q_p = Q - quantities[x_i]
            if q_p > 0:
                arg_min_1 = min(g[x_i], key=g[x_i].get)
                p[l][x_i] = arg_min_1
                f[l][x_i] = g[x_i][arg_min_1]
                del g[x_i][arg_min_1]
                arg_min_2 = min(g[x_i], key=g[x_i].get)
                phi[l][x_i] = g[x_i][arg_min_2]
                q_route[l][x_i] = g_type[x_i][arg_min_1]
                q_route_2[l][x_i] = g_type[x_i][arg_min_2]
    return (f, phi, p, q_route, q_route_2)

def construct_q_paths_debug(h,truck_capacity,N,distance,values,values_pos,quantities,direction):

    # Initialize the routes
    f = {}; phi = {}; p = {}; q_route = {}; q_route_2 = {};
    for l in range(len(values)):
        f[l] = {}; phi[l] = {}; p[l] = {}; q_route[l] = {}; q_route_2[l] = {};
        for n in N:
            f[l][n] = float('inf'); phi[l][n] = float('inf'); p[l][n] = float('inf');
            q_route[l][n] = []; q_route_2[l][n] = [];

    # Initialize the routes
    for n in N:
        q = quantities[n]
        if q<= truck_capacity:
            l = values_pos[q]
            if direction == 'left':
                f[l][n] = distance[h][n]
            else:
                f[l][n] = distance[n][h]
            p[l][n] = h
            q_route[l][n] = [h, n]

    # Calculate the recursion
    for l in range(len(values)):
        Q = values[l]
        g = {}
        g_type = {}
        for x_i in N:
            q_p = Q - quantities[x_i]
            g[x_i] = {}
            g_type[x_i] = {}
            if q_p>0:
                l_p = values_pos[q_p]
                for x_j in N:
                    if x_i!=x_j:
                        if p[l_p][x_j]!=x_i:
                            g[x_i][x_j] = f[l_p][x_j] + distance[x_j][x_i]
                            g_type[x_i][x_j] = (0,l_p)
                        else:
                            g[x_i][x_j] = phi[l_p][x_j] + distance[x_j][x_i]
                            g_type[x_i][x_j] = (1,l_p)
        for x_i in N:
            q_p = Q - quantities[x_i]
            if q_p > 0:
                min_1 = float('inf')
                min_2 = float('inf')
                arg_min_1 = float('inf')
                arg_min_2 = float('inf')
                for x_j in N:
                    if (x_i!=x_j):
                        value = g[x_i][x_j];
                        if (value<min_1):
                            min_2 = min_1
                            min_1 = value
                            arg_min_2 = arg_min_1
                            arg_min_1 = x_j
                        elif (value<min_2):
                            min_2 = value
                            arg_min_2 = x_j

                p[l][x_i] = arg_min_1
                f[l][x_i] = g[x_i][arg_min_1]
                phi[l][x_i] = g[x_i][arg_min_2]
                coord = g_type[x_i][arg_min_1]
                coord_2 = g_type[x_i][arg_min_2]
                q_route[l][x_i] = (q_route[coord[1]][arg_min_1] + [x_i]) if (coord[0] == 0) else (q_route_2[coord[1]][arg_min_1] + [x_i])
                q_route_2[l][x_i] = (q_route[coord_2[1]][arg_min_2] + [x_i]) if (coord_2[0] == 0) else (q_route_2[coord_2[1]][arg_min_2] + [x_i])
    return (f, phi, p, q_route, q_route_2)



# Combines paths to construct routes
def construct_q_routes(h,truck_capacity,N,distance,values,values_pos,quantities):
    f_l, phi_l, p_l, q_route_l, q_route_2_l = construct_q_paths(h,truck_capacity,N,distance,values,values_pos,quantities,'left')
    f_r, phi_r, p_r, q_route_r, q_route_2_r = construct_q_paths(h,truck_capacity,N,distance,values,values_pos,quantities,'right')

    # Initialize the routes
    psi = {}; psi_route = {};
    for l in range(len(values)):
        psi[l] = {}; psi_route[l] = {};
        for n in N:
            psi[l][n] = {}; psi_route[l][n] = {};

    for l in range(len(values)):
        for n in N:
            min_w = quantities[n]
            max_w = (values[l] + quantities[n])
            min_route = []
            min_val = float('inf')
            for l_1 in range(len(values)):
                q = values[l_1]
                if q>=min_w and q<max_w:
                    l_2 = values_pos[values[l] + quantities[n]-q]
                    if p_l[l_1][n]!=p_r[l_2][n] or (p_l[l_1][n] == p_r[l_2][n] and p_l[l_1][n]==h):
                        val = f_l[l_1][n]+f_r[l_2][n]
                        route = q_route_l[l_1][n]+list(reversed(q_route_r[l_2][n]))[1:]
                    else:
                        if f_l[l_1][n]+phi_r[l_2][n]<phi_l[l_1][n]+f_r[l_2][n]:
                            val = f_l[l_1][n]+phi_r[l_2][n]
                            route = q_route_l[l_1][n]+list(reversed(q_route_2_r[l_2][n]))[1:]
                        else:
                            val = phi_l[l_1][n]+f_r[l_2][n]
                            route = q_route_2_l[l_1][n]+list(reversed(q_route_r[l_2][n]))[1:]
                    if val<min_val:
                        min_val = val
                        min_route = route
            psi[l][n] = min_val
            psi_route[l][n] = min_route
    return (psi, psi_route)

# Combines paths to construct routes
def construct_q_routes_debug(h,truck_capacity,N,distance,values,values_pos,quantities):
    f_l, phi_l, p_l, q_route_l, q_route_2_l = construct_q_paths_debug(h,truck_capacity,N,distance,values,values_pos,quantities,'left')
    f_r, phi_r, p_r, q_route_r, q_route_2_r = construct_q_paths_debug(h,truck_capacity,N,distance,values,values_pos,quantities,'right')

    # Initialize the routes
    psi = {}; psi_route = {};
    for l in range(len(values)):
        psi[l] = {}; psi_route[l] = {};
        for n in N:
            psi[l][n] = {}; psi_route[l][n] = {};

    for l in range(len(values)):
        for n in N:
            min_w = quantities[n]
            max_w = (values[l] + quantities[n])
            min_route = []
            min_val = float('inf')
            min_coord = ()
            for l_1 in range(len(values)):
                q = values[l_1]
                if q>=min_w and q<max_w:
                    l_2 = values_pos[values[l] + quantities[n]-q]
                    if p_l[l_1][n]!=p_r[l_2][n] or (p_l[l_1][n] == p_r[l_2][n] and p_l[l_1][n]==h):
                        val = f_l[l_1][n]+f_r[l_2][n]
                        coord = (0,l_1,l_2)
                    else:
                        if f_l[l_1][n]+phi_r[l_2][n]<phi_l[l_1][n]+f_r[l_2][n]:
                            val = f_l[l_1][n]+phi_r[l_2][n]
                            coord = (1,l_1,l_2)
                        else:
                            val = phi_l[l_1][n]+f_r[l_2][n]
                            coord = (2,l_1,l_2)
                    if val<min_val:
                        min_val = val
                        min_coord = coord
            psi[l][n] = min_val
            if len(min_coord)>0:
                if min_coord[0] == 0:
                    min_route = q_route_l[min_coord[1]][n]+list(reversed(q_route_r[min_coord[2]][n]))[1:]
                elif min_coord[0] == 1:
                    min_route = q_route_l[min_coord[1]][n]+list(reversed(q_route_2_r[min_coord[2]][n]))[1:]
                elif min_coord[0] == 2:
                    min_route = q_route_2_l[min_coord[1]][n]+list(reversed(q_route_r[min_coord[2]][n]))[1:]
            psi_route[l][n] = min_route
    return (psi, psi_route)



# Calculate lower bound of the optimization problem, along with auxiliary quantities
def lower_bound(H,capacities,N,quantities,distance,mu,lamb):
    # Calculate the bound of the problem and the routes that help us reach that bound
    psi = {};
    psi_route = {};
    b = {};
    W = {}
    for h in H:
        truck_capacity = capacities[h]
        values = possible_values(quantities,capacities[h])
        W[h] = values
        values_pos = dict(zip(values,list(range(len(values)))))
        psi_h, psi_route_h = construct_q_routes_debug(h,truck_capacity,N,distance,values,values_pos,quantities)
        psi[h] = psi_h
        psi_route[h] = psi_route_h
        b[h] = {}
        for l in range(len(values)):
            b[h][l] = {}
            for n in N:
                b[h][l][n] = (psi[h][l][n]-mu[h])*quantities[n]/values[l]

    min_b = {}
    arg_l = {}
    for h in H:
        min_b[h] = {}
        arg_l[h] = {}
        for n in N:
            min_b[h][n] = np.min(np.array([b[h][l][n] for l in range(len(W[h]))]))
            arg_l[h][n] = np.argmin(np.array([b[h][l][n] for l in range(len(W[h]))]))

    min_min_b = {}
    arg_h = {}
    arg_arg_l = {}
    arg_route = {}
    for n in N:
        i = True
        min_min_b[n] = np.min(np.array([min_b[h][n] for h in H]))
        arg_h[n] = H[np.argmin(np.array([min_b[h][n] for h in H]))]
        arg_arg_l[n] = arg_l[arg_h[n]][n]
        arg_route[n] = psi_route[arg_h[n]][arg_l[arg_h[n]][n]][n]

    # Compute theta and rho
    # First we calculate the number of times each node is visited
    visits = {}
    for n in N:
        visits[n] = {}
        for n_2 in N:
            visits[n][n_2] = np.sum([n_2 == i for i in arg_route[n]])

    theta = {}
    for j in N:
        theta[j] = 0
        for i in N:
            theta[j] = theta[j] + (quantities[i]+0.0)/W[arg_h[i]][arg_arg_l[i]]*visits[i][j]

    rho = {}
    for h in H:
        rho[h] = 0
    for i in N:
        rho[arg_h[i]] = rho[arg_h[i]] + (quantities[i]+0.0)/W[arg_h[i]][arg_arg_l[i]]

    # Construct the duality vectors
    u = {}
    for n in N:
        u[n] = min_min_b[n] + lamb[n]

    z_lb = np.sum([u[n] for n in N]) + np.sum([mu[h] for h in H])

    return (z_lb, theta, rho, u)



def optimize_lower_bound(iterations, z_ub, epsilon, H,capacities,N,quantities,distance_dictionary):

    # Initialize the parameters
    lamb = {}
    mu = {}
    for j in N:
        lamb[j] = 0
    for h in H:
        mu[h] = 0
    max_val = -float('inf')
    values = []

    for i in range(iterations):
        distance = reduced_cost_dict(lamb, distance_dictionary, N)

        z_lb, theta, rho, u = lower_bound(H,capacities,N,quantities,distance,mu,lamb)
        values.append(z_lb)
        if z_lb > max_val:
            max_val = copy.deepcopy(z_lb)
            u_opt = copy.deepcopy(u)
            v_opt = copy.deepcopy(mu)
            lamb_opt = copy.deepcopy(lamb)

        # Compute the new parameters
        gamma = (z_ub - z_lb)/(np.sum((np.array(theta.values())-1)**2) + np.sum((np.array(rho.values())-1)**2))
        print(z_lb)
        print(lamb)
        print(mu)
        # New lambda
        for j in N:
            lamb[j] = lamb[j] - epsilon*gamma*(theta[j]-1)
        for h in H:
            mu[h] = np.min([mu[h] - epsilon*gamma*(rho[h]-1),0])


        if np.sum(np.abs(lamb.values())) > 10**16:
            raise ValueError('Lambda exploding')

        # Rule for updating epsilon
        if len(values)>=7:
            grad = [values[i+1]-values[i] for i in range(len(values)-1)]
            jumps = [np.sign(grad[i+1])!=np.sign(grad[i]) for i in range(len(grad)-1)]
            if np.sum(jumps[len(jumps)-5:len(jumps)])>=3:
                epsilon = epsilon/1.5
                print ('new epsilon:%f' % epsilon)
                values = []
            if np.sum(np.array(grad[len(grad)-5:len(grad)])>0) >= 5:

                epsilon = epsilon*1.2
                print ('new epsilon:%f' % epsilon)
                values = []
        if (np.array_equal(np.array(theta.values()), np.ones(len(N))) and np.array_equal(np.array(rho.values()), np.ones(len(H)))):
            print("reached zero gradient")
            return (max_val,u_opt,v_opt,lamb_opt)
    return (max_val,u_opt,v_opt,lamb_opt)

def optimize_lower_bound_c(iterations, z_ub, epsilon, H,capacities,N,quantities,distance_mat):
    H_ = (np.array(range(len(H)))+len(N)).astype(int)
    N_ = np.array(range(len(N))).astype(int)
    capacities_ = np.array([capacities[h] for h in H]).astype(int)
    quantities_ = np.array([quantities[n] for n in N]).astype(int)
    geo_distance_ = distance_mat.astype("float64")
    result = cpp_lower_bounds.lower_bound_optimizer_M1_(iterations, z_ub, epsilon, H_, capacities_, N_, quantities_, geo_distance_)

    max_val = result["z_lb"]
    u_opt = result["u"]
    v_opt = result["v"]
    lamb_opt = result["lamb"]

    return (max_val,u_opt,v_opt,lamb_opt)

# We solve the GENPATH problem
def GENPATH(Delta, gamma, h, capacity, N, quantities, distance, direction):
    P = {}
    T = {}
    for k in N + [h]:
        P[k] = []
        T[k] = []
    T[h] = [{'path':[h],'cost':0,'lower_bound':0, "load":0, 'end':h}]

    count_paths = 0
    while True:
        costs = {}
        for k in T.keys():
            if len(T[k])>0:
                costs[k] = T[k][0]['cost']
        if len(costs) == 0:
            break
        min_costs = min(costs, key = costs.get)

        p_star = T[min_costs].pop(0)
        if not min_costs in P.keys():
            P[min_costs] = []
        P[min_costs].append(p_star)
        count_paths += 1
        # If too many paths, stop
        if count_paths==Delta:
            break
        # If path violates capacity, go to the next one
        if p_star['load'] > capacity/2.0:
            continue
        for n in N:
            if not (n in p_star['path']):
                if direction == 'right':
                    new_p = {'path':p_star['path'] + [n], 'cost':p_star['cost'] + distance[n][p_star['end']],
                            'lower_bound':p_star['cost'] + distance[n][p_star['end']], 'load': p_star['load'] + quantities[n],
                            'end':n}
                elif direction == 'left':
                    new_p = {'path':p_star['path'] + [n], 'cost':p_star['cost'] + distance[p_star['end']][n],
                            'lower_bound':p_star['cost'] + distance[p_star['end']][n], 'load': p_star['load'] + quantities[n],
                            'end':n}
                # Check if the new path has a cost too high
                if new_p['lower_bound'] >= gamma:
                    continue
                # Check if the new path has a load too high
                if new_p['load'] > capacity:
                    continue
                # Check if this new path is dominated by any path in P
                dominated = False
                for p in P[n]:
                    if (p['end'] == new_p['end']) and (p['cost'] <= new_p['cost']) and (set(p['path']) == set(new_p['path'])):
                        dominated = True
                        break
                if dominated:
                    continue
                # Check if the path is dominated by any path in T
                insertion_index = 0
                for i,p in enumerate(T[n]):
                    if (p['end'] == new_p['end']) and (p['cost'] <= new_p['cost']) and (set(p['path']) == set(new_p['path'])):
                        dominated = True
                        break
                    if (p['cost'] > new_p['cost']):
                        break
                    insertion_index = i+1
                if dominated:
                    continue
                # Append the path
                T[n].insert(insertion_index, new_p)
                # Delete dominated elements
                j = insertion_index + 1
                while j<len(T[n]):
                    p = T[n][j]
                    if (p['end'] == new_p['end']) and (p['cost'] > new_p['cost']) and (set(p['path']) == set(new_p['path'])):
                        T[n].pop(j)
                    else:
                        j += 1
    return P

def GENROUTE(Delta, gamma, h, capacity, N, quantities, distance):
    P_l = GENPATH(Delta, gamma, h, capacity, N, quantities, distance, direction = 'left')
    P_r = GENPATH(Delta, gamma, h, capacity, N, quantities, distance, direction = 'right')

    T = {}
    R = {}
    added = {}
    for n in N:
        added[n] = set((-1,-1))
        if len(P_l[n])>0 and len(P_r[n])>0:
            T[n] = [[(0,0),P_l[n][0]['cost']+P_r[n][0]['cost']]]
            added[n].add((0,0))
        else:
            T[n] = []
        R[n] = []

    valid_v = [0,0,0,0]
    while True:
        # ipdb.set_trace()
        # Calculate costs
        costs = {}
        for n in N:
            if len(T[n])>0:
                costs[n] = T[n][0][1]
        if len(costs) == 0:
            break
        min_costs_n = min(costs, key = costs.get)
        min_cost = costs[min_costs_n]
        indices = T[min_costs_n].pop(0)[0]
        path_l = P_l[min_costs_n][indices[0]]
        path_r = P_r[min_costs_n][indices[1]]
        if min_cost> gamma:
            break
        total_load = path_l['load'] + path_r['load'] - quantities[min_costs_n]
        valid = True
        if total_load > capacity:
            valid = False
            valid_v[0] = valid_v[0]+1

        elif (np.min([path_l['load'],path_r['load']]) < total_load/2.0 or
            np.max([path_l['load'],path_r['load']]) > total_load/2.0+quantities[min_costs_n]):
            valid = False
            valid_v[1] = valid_v[1]+1

        elif (set(path_l['path']).intersection(set(path_r['path'])) != set([h,min_costs_n])):
            valid = False
            valid_v[2] = valid_v[2]+1
        else:
            for n in N:
                for r in R[n]:
                    if set(r['path']) == set(path_l['path']+path_r['path']):
                        valid = False
                        valid_v[3] = valid_v[3] + 1
                        break
                if not valid:
                    break
        if valid:
            R[min_costs_n].append({'path':path_l['path'][0:(len(path_l['path'])-1)]+list(reversed(path_r['path'])),
                                  'cost':path_l['cost']+path_r['cost'],
                                  'load':total_load,
                                  'median':min_costs_n,
                                  'indices':indices})

        new_route_1 = (indices[0]+1,indices[1])
        new_route_2 = (indices[0],indices[1]+1)
        # If routes do not exist, transform them into the first route
        if (indices[0]+1 >= len (P_l[min_costs_n])):
            new_route_1 = (0,0)
        if (indices[1]+1 >= len (P_r[min_costs_n])):
            new_route_2 = (0,0)
        new_routes = [new_route_1,new_route_2]
        new_costs = [P_l[min_costs_n][new_routes[0][0]]['cost']+P_r[min_costs_n][new_routes[0][1]]['cost'],
                     P_l[min_costs_n][new_routes[1][0]]['cost']+P_r[min_costs_n][new_routes[1][1]]['cost']]
        min_cost = np.min(new_costs)
        max_cost = np.max(new_costs)
        min_route = new_routes[np.argmin(new_costs)]
        max_route = new_routes[(np.argmin(new_costs)+1)%2]
        insert_index = 0
        # Check if the route has been added previously
        if not min_route in added[min_costs_n]:
            for i in range(len(T[min_costs_n])):
                cost = T[min_costs_n][i][1]
                if min_cost<cost:
                    break
                insert_index+=1
            T[min_costs_n].insert(insert_index,[min_route,min_cost])
            insert_index +=1
            added[min_costs_n].add(min_route)
        # Check if the route has been added previously
        if not max_route in added[min_costs_n]:
            for i in range(insert_index, len(T[min_costs_n])):
                cost = T[min_costs_n][i][1]
                if max_cost<cost:
                    break
                insert_index+=1
            T[min_costs_n].insert(insert_index,[max_route,max_cost])
            added[min_costs_n].add(max_route)

    # Verify that the routes are not always empty
    for n in N:
        if len(R[n]) == 0:
            R[n] = [{'path':[h,n,h], 'cost':distance[h][n] + distance[n][h], 'load': quantities[n]}]
    return R


def GENPATH_c(Delta, gamma, h, capacity, N, quantities, distance_mat, direction):
    N_ = np.array(range(len(N))).astype(int)
    quantities_ = np.array([quantities[n] for n in N]).astype(int)

    results_c = cpp_lower_bounds.GENPATH_(Delta, gamma, h, capacity, N_, quantities_, distance_mat, direction)
    # We make sure to print in the same format
    results = {}
    for i,n in enumerate(N):
        results[n] = results_c[i]
        for path in results[n]:
            del path['nodes']
            path['end'] = 'n_'+str(path['end']) if (path['end']<len(N)) else 'h_'+str(path['end']-len(N))
            for j,f in enumerate(path['path']):
                path['path'][j] = 'n_'+str(f) if (f<len(N)) else 'h_'+str(f-len(N))
    results['h_'+str(h-len(N))] = results_c[len(N)]
    for path in results['h_'+str(h-len(N))]:
        del path['nodes']
        path['end'] = 'n_'+str(path['end']) if (path['end']<len(N)) else 'h_'+str(path['end']-len(N))
        for j,f in enumerate(path['path']):
            path['path'][j] = 'n_'+str(f) if (f<len(N)) else 'h_'+str(f-len(N))

    return results

def GENROUTE_c(Delta, gamma, h, capacity, N, quantities, distance_mat, geo_distance):
    N_ = np.array(range(len(N))).astype(int)
    quantities_ = np.array([quantities[n] for n in N]).astype(int)

    results_c = cpp_lower_bounds.GENROUTE_(Delta, gamma, h, capacity, N_, quantities_, distance_mat, geo_distance)
    # We make sure to print in the same format
    results = results_c
    for path in results:
        #path['median'] = 'n_'+str(path['median']) if (path['median']<len(N)) else 'h_'+str(path['median']-len(N))
        #path['indices'] = (path['index_l'],path['index_r'])
        del path['median']
        del path['index_r']
        del path['index_l']
        #del path['truck']
        #del path['geo_cost']

        #for j,f in enumerate(path['path']):
        #    path['path'][j] = 'n_'+str(f) if (f<len(N)) else 'h_'+str(f-len(N))

    return results

def lower_bound_M2(Routes, H,capacities,N,quantities,distance,mu,lamb,geo_distance):

    # Given the routes, calculate the dual variables
    b = {}
    b_route = {}
    for h in H:
        b[h] = {}
        b_route[h] = {}
        for n in N:
            b[h][n] = float('inf')
            b_route[h][n] = []

    for h in H:
        for n in Routes[h].keys():
            for r in Routes[h][n]:
                nodes = r['path'][1:(len(r['path'])-1)]
                cost_r = 0
                for i in range(len(r['path'])-1):
                    cost_r += distance[r['path'][i]][r['path'][i+1]]
                for node in nodes:
                    new_b = (cost_r-mu[h])/r['load']*quantities[node]
                    if new_b<b[h][node]:
                        b[h][node] = new_b
                        b_route[h][node] = r

    min_b = {}
    min_b_route = {}
    for n in N:
        min_b[n] = float('inf')
        min_b_route[n] = {}
        for h in H:
            if b[h][n]<min_b[n]:
                min_b[n] = b[h][n]
                min_b_route[n] = {'path':b_route[h][n]['path'], 'truck':h, 'load': b_route[h][n]['load']}

    # Construct the update vectors
    rho = {}
    for h in H:
        rho[h] = 0
    for n in N:
        route = min_b_route[n]
        rho[route['truck']] += (quantities[n]+0.0)/route['load']

    theta = {}
    for n in N:
        theta[n] = 0
    for n in N:
        route = min_b_route[n]
        for v in route['path']:
            if v in N:
                theta[v] += (quantities[n]+0.0)/route['load']

    # Construct the duality vectors
    u = {}
    for n in N:
        u[n] = min_b[n] + lamb[n]

    z_lb = np.sum([u[n] for n in N]) + np.sum([mu[h] for h in H])

    # Update the costs of the routes
    for h in H:
        for n in Routes[h].keys():
            for r in Routes[h][n]:
                nodes = r['path'][1:(len(r['path'])-1)]
                geo_cost = 0
                u_cost = 0
                for i in range(len(r['path'])-1):
                    geo_cost += geo_distance[r['path'][i]][r['path'][i+1]]
                for node in nodes:
                    u_cost += u[node]
                r['reduced_cost'] = geo_cost - u_cost - mu[h]
                r['geo_cost'] = geo_cost

    return (z_lb, theta, rho, u)

def maximize_lower_bound_M2(Routes, lamb, mu, geo_distance, z_ub, iterations, epsilon,N,H, capacities, quantities):
    max_val = -float('inf')
    values = []

    for i in range(iterations):
        distance = reduced_cost_dict(lamb, geo_distance,N)

        z_lb, theta, rho, u = lower_bound_M2(Routes, H,capacities,N,quantities,distance,mu,lamb,geo_distance)
        values.append(z_lb)
        if z_lb > max_val:
            max_val = copy.deepcopy(z_lb)
            u_opt = copy.deepcopy(u)
            v_opt = copy.deepcopy(mu)
            lamb_opt = copy.deepcopy(lamb)

        # Compute the new parameters
        gamma = (z_ub - z_lb)/(np.sum((np.array(theta.values())-1)**2) + np.sum((np.array(rho.values())-1)**2))
        # New lambda
        for j in N:
            lamb[j] = lamb[j] - epsilon*gamma*(theta[j]-1)
        for h in H:
            mu[h] = np.min([mu[h] - epsilon*gamma*(rho[h]-1),0])
        print(z_lb)

        if np.sum(np.abs(lamb.values())) > 10**16:
            raise ValueError('Lambda exploding')

        # Rule for updating epsilon
        if len(values)>=7:
            grad = [values[i+1]-values[i] for i in range(len(values)-1)]
            jumps = [np.sign(grad[i+1])!=np.sign(grad[i]) for i in range(len(grad)-1)]
            if np.sum(jumps[len(jumps)-5:len(jumps)])>=3:
                epsilon = epsilon/1.5
                print ('new epsilon:%f' % epsilon)
                values = []
            if np.sum(np.array(grad[len(grad)-5:len(grad)])>0) >= 5:

                epsilon = epsilon*1.2
                print ('new epsilon:%f' % epsilon)
                values = []
        if (np.array_equal(np.array(theta.values()), np.ones(len(N))) and np.array_equal(np.array(rho.values()), np.ones(len(H)))):
            print("reached zero gradient")
            return (max_val,u_opt,v_opt,lamb_opt)

    return (max_val,u_opt,v_opt,lamb_opt)

def lower_bound_optimizer_M2_c(sub_iterations, z_ub, Delta, Delta_zero, gamma, gamma_zero, epsilon, H, capacities, N, quantities, geo_distance, mu, lamb):

    H_ = (np.array(range(len(H)))+len(N)).astype(int)
    N_ = np.array(range(len(N))).astype(int)
    capacities_ = np.array([capacities[h] for h in H]).astype(int)
    quantities_ = np.array([quantities[n] for n in N]).astype(int)
    geo_distance_ = geo_distance.astype("float64")
    mu_ = np.array([mu[h] for h in H])
    lamb_ = np.array([lamb[n] for n in N])

    result = cpp_lower_bounds.optimize_lower_bound_M2_(sub_iterations, z_ub, Delta, Delta_zero, gamma, gamma_zero, epsilon, H_, capacities_, N_, quantities_, geo_distance_, mu_, lamb_)

    max_val = result["z_lb"]
    u_opt = result["u"]
    v_opt = result["v"]
    lamb_opt = result["lamb"]

    return (max_val,u_opt,v_opt,lamb_opt)

def lower_bound_optimizer_M2(sub_iterations, z_ub, Delta, Delta_zero, gamma, gamma_zero, epsilon, H, capacities, N, quantities, geo_distance, mu, lamb):

    reduced_distance = reduced_cost_dict_complete(lamb, mu, geo_distance,N,H)
    # Calculate the routes of minimum cost
    Routes = {}
    for h in H:
        Routes[h] = GENROUTE(Delta, gamma, h, capacities[h], N, quantities, reduced_distance)

    max_val,u_opt,v_opt,lamb_opt = maximize_lower_bound_M2(Routes, lamb, mu, geo_distance, z_ub, sub_iterations, epsilon, N,H, capacities, quantities)

    return (max_val,u_opt,v_opt,lamb_opt)

def construct_lower_bound_c(iterations_grad_m1,iterations_grad_m2,iterations_m2,z_ub,Delta,Delta_zero,Delta_final,gamma,gamma_zero,gamma_final,epsilon,H,capacities,N,quantities,geo_distance,n_trucks):
    H_ = (np.array(range(len(H)))+len(N)).astype(int)
    N_ = np.array(range(len(N))).astype(int)
    capacities_ = np.array([capacities[h] for h in H]).astype(int)
    quantities_ = np.array([quantities[n] for n in N]).astype(int)
    geo_distance_ = geo_distance.astype("float64")
    n_trucks_ = np.array([n_trucks[h] for h in H]).astype(int)

    result = cpp_lower_bounds.construct_lower_bound_(iterations_grad_m1,iterations_grad_m2,iterations_m2,z_ub,Delta,Delta_zero,Delta_final,gamma,gamma_zero,gamma_final,epsilon,H_,capacities_,N_,quantities_,geo_distance_,n_trucks_)


    return (result)

def primal_solver(solution, len_N, H, quantities, capacities, n_trucks, time):
    print("The lower bound is %f"%solution["z_lb"])
    reduced_routes =[r for truck_routes in solution["routes"] for r in truck_routes]
    print("There are %d total routes"%len(reduced_routes))
    truck_routes={}
    farmer_routes={}
 
    for farmer in range(len_N):
        farmer_routes[farmer] = [i for i,route in enumerate(reduced_routes) if farmer in route["path"]]
        print(farmer)
    for truck in range(len_N, len_N + len(H)):
        truck_routes[truck] = [i for i,route in enumerate(reduced_routes) if truck==route["truck"]]
        print(truck)
    model = Model()
    # Add variables to modes
    variables = []
    for i,route in enumerate(reduced_routes):
        variables.append(model.addVar(obj=route['geo_cost'], vtype=GRB.BINARY,
                                      name='route_'+str(i)))
    model.update()
    # Add farmer constraints
    for farmer, routes in farmer_routes.iteritems():
        model.addConstr(quicksum([variables[k] for k in routes])>=1)
    ## Add route constraints
    for truck, routes in truck_routes.iteritems():
        model.addConstr(quicksum([variables[k] for k in routes])<=n_trucks['h_'+str(truck-len_N)])
    model.update()
    model.setParam('TimeLimit', time)
    model.update()
 
    model.optimize()

    if True:
        ## Extract the optimal Routes
        routes_chosen = []
        for v in model.getVars():
            name = v.VarName
            value = v.x
            if value == 1.0:
                r_number = int(name.split('_')[1])
                t_number = int(name.split('_')[2])
                route_chosen = copy.deepcopy(reduced_routes[r_number])
                route_chosen["truck"] = t_number
                routes_chosen.append(route_chosen)

        # Format them
        routes_formatted = {h:[] for h in H}
        for r in routes_chosen:
            h = 'h_'+str(r["truck"] - len_N)
            new_route = {}
            new_route["route"] = ['n_'+str(i) for i in r['path'][1:(len(r['path'])-1)]]
            new_route["load"] = np.sum([quantities[n] for n in new_route["route"]])
            new_route["max_load"] = capacities[h]
            new_route["index"] = r["truck"]
            routes_formatted[h].append(new_route)
    return routes_formatted

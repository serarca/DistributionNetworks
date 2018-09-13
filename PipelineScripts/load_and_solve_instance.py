# This script takes the data for each instance of our temporal and
# spatial problem and makes it into pickle files
import datetime
from dateutil import parser
import math
from dateutil.parser import parse
import pandas as pd
import numpy as np
from shapely.geometry import Point
import time
import sys

sys.path.insert(0, 'Code/VRPEngine/pyCode')
sys.path.insert(0, 'Code/VRPEngine/C++Engine')

instance = sys.argv[1]

import solver as solver
import VRPClass
import pickle
import lower_bound

# This script loads a given instance and calculates its optimal solution
# calling the C++ library

# Unpickle
import pickle
file_dict = pickle.load( open( instance + ".p", "rb" ) )
print(file_dict)
locals().update(file_dict)
# Get data
print("There is a total of %d trucks"% len(H))
print("There is a total of %d farmers"% len(N))

# Let's transform quantities in integers
quantities = {k:int(10*quantities[k]) for k in quantities.keys()}
capacities = {k:int(10*capacities[k]) for k in capacities.keys()}

print(quantities)
print(capacities)
geo_distance = distance

# Parameters of the algorithms
Delta = 3000
gamma = 20000
iterations_m2 = 2
z_ub = 500000


# Parameters for the search algorithm that will not change
m = 1
iterations_grad_m1 = 200
iterations_grad_m2 = 100
Delta_zero = 1000
Delta_final = Delta
gamma_zero = -10**(-14) * z_ub
gamma_final = gamma
epsilon = 0.1
time_limit = 200


result = lower_bound.construct_lower_bound_c(
    iterations_grad_m1,iterations_grad_m2,
    iterations_m2,z_ub,Delta,Delta_zero,
    Delta_final,gamma,gamma_zero,gamma_final,
    epsilon,H,capacities,N,quantities,geo_distance,n_trucks)

routes = lower_bound.primal_solver(result[iterations_m2],len(N),H, quantities, capacities, n_trucks, time_limit)

vrp = VRPClass.VRP(H, N, H_p, N_p, quantities, capacities, type_dist, M = M, M_p = M_p, distance_matrix = distance)
z_ub = np.sum([vrp.distance_path([h]+route['route']+[h]) for h in H for route in routes[h]])
z_lb = result[iterations_m2]['z_lb']
z_gap = (z_ub - z_lb)/z_ub

result_dict = {"routes":routes,
              "mapping":mapping,
              "z_ub":z_ub,
              "z_lb":z_lb,
              "z_gap":z_gap}
import json
with open( instance + "_result.json", 'wb') as outfile:
    json.dump(result_dict, outfile)

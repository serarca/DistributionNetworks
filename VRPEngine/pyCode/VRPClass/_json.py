import json
import numpy as np
from ._primal_solution import *



@classmethod
def from_json(cls, file_name):
	with open(file_name) as json_file:  
	    data = json.load(json_file)
	    H = data['H']
	    N = data['N']
	    H_p = data['H_p']
	    N_p = data['N_p']
	    capacities = data['capacities']
	    quantities = data['quantities']
	    type_dist = 'road'
	    distance_matrix = np.array(data['distances'])
	    mapping = data['mapping']
	    n_trucks = data['n_trucks']
	    penalties = data['penalties']


	vrp = cls(H, N, H_p, N_p, quantities, capacities, type_dist, distance_matrix = distance_matrix)
	vrp.mapping = mapping
	vrp.n_trucks = n_trucks
	vrp.penalties = penalties

	return vrp

def primal_from_json(file_name):
	with open(file_name) as json_file: 
	    data = json.load(json_file)
	    

	primal = PrimalSolution(data["routes"], data["z_lb"])

	return primal


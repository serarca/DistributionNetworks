import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import math
import random
from gurobipy import *
import geopy.distance
import pandas as pd
from tsp_nn import NearestNeighbourSolver
from tsp_ni import NearestInsertionSolver
import copy
from scipy.spatial import distance as distance_lib
import lower_bound
import road_distance



# This class holds a VRP problem
class VRP:
	def __init__(self, H, N, H_p, N_p, quantities, capacities, type_dist, M = False, M_p = False, distance_matrix = None):
		self.V = N+H
		self.H = H
		self.N = N
		self.H_p = H_p
		self.N_p = N_p
		self.quantities = quantities
		self.capacities = capacities
		self.type_dist = type_dist
		self.M_p = M_p
		self.M = M
		self.time = None

		# Find closest mill to each node
		if M:
			self.closest_mill = {}
			for i in range(len(self.N)):
				try:
					self.closest_mill[N[i]] = M[distance_lib.cdist([N_p[i,:]], M_p).argmin()]
				except:
					ipdb.set_trace()
		# Create distance matrices
		if type(distance_matrix) is np.ndarray:
			self.distance = distance_matrix
		else:
			if type_dist == 'geo':
				self.geo_distance = self.geo_distance_matrix(pd.DataFrame(data = np.concatenate([N_p, H_p]), index = np.concatenate([self.N,self.H]), columns = ['lat','lon']))
				self.distance = self.geo_distance
			if type_dist == 'euclid':
				self.euclid_distance = self.euclid_distance_matrix(pd.DataFrame(data = np.concatenate([N_p, H_p]), index = np.concatenate([self.N,self.H]), columns = ['lat','lon']))
				self.distance = self.euclid_distance
			if type_dist == 'manhattan':
				self.manhattan_distance = self.manhattan_distance_matrix(pd.DataFrame(data = np.concatenate([N_p, H_p]), index = np.concatenate([self.N,self.H]), columns = ['lat','lon']))
				self.distance = self.manhattan_distance
			if type_dist == 'road':
				aux = self.road_distance_matrix(pd.DataFrame(data = np.concatenate([N_p, H_p]), index = np.concatenate([self.N,self.H]), columns = ['lat','lon']))
				self.road_distance = aux[0]
				self.road_time = aux[1]
				self.distance = self.road_distance
				self.time = self.road_time
		self.positions = dict(zip(np.concatenate([N,H]),list(range(len(N+H)))))
		# Create distance dictionary
		distance_dict = {}
		for v1 in self.V:
			distance_dict[v1] = {}
			for v2 in self.V:
				distance_dict[v1][v2] = self.distance[self.positions[v1],self.positions[v2]]
		self.distance_dict = distance_dict
		if (self.time is not None):
			time_dict = {}
			for v1 in self.V:
				time_dict[v1] = {}
				for v2 in self.V:
					time_dict[v1][v2] = self.time[self.positions[v1],self.positions[v2]]
			self.time_dict = time_dict


	# Construct geographical distances matrix
	def geo_distance_matrix(self,df,mills = False):
		df = df.reset_index(drop=True).copy()
		d = np.zeros([len(df),len(df)])
		for i,row_i in df.iterrows():
			for j,row_j in df.iterrows():
				if i!=j:
					d[i,j] = geopy.distance.vincenty(df.iloc[i].values,df.iloc[j].values).m
					if (self.M):
						if ((self.V[j] in self.H) and (self.V[i] in self.N)):
							d[i,j] = (geopy.distance.vincenty(df.iloc[i].values,self.get_entry(*str.split(self.closest_mill[self.V[i]],'_'))).m +
									  geopy.distance.vincenty(self.get_entry(*str.split(self.closest_mill[self.V[i]],'_')),df.iloc[j].values).m)
		return d

	# Construct euclid distances matrix
	def euclid_distance_matrix(self,df):
		df = df.reset_index(drop=True).copy()
		d = np.zeros([len(df),len(df)])
		for i,row_i in df.iterrows():
			for j,row_j in df.iterrows():
				if i!=j:
					d[i,j] = np.sqrt(np.sum((df.iloc[i].values-df.iloc[j].values)**2))
					if (self.M):
						if ((self.V[j] in self.H) and (self.V[i] in self.N)):
							d[i,j] = (np.sqrt(np.sum((df.iloc[i].values-self.get_entry(*str.split(self.closest_mill[self.V[i]],'_')))**2)) +
									  np.sqrt(np.sum((self.get_entry(*str.split(self.closest_mill[self.V[i]],'_'))-df.iloc[j].values)**2)))

		return d

	# Construct manhattan distances matrix
	def manhattan_distance_matrix(self,df):
		df = df.reset_index(drop=True).copy()
		d = np.zeros([len(df),len(df)])
		for i,row_i in df.iterrows():
			for j,row_j in df.iterrows():
				if i!=j:
					d[i,j] = np.sum(np.abs(df.iloc[i].values-df.iloc[j].values))
					if (self.M):
						if ((self.V[j] in self.H) and (self.V[i] in self.N)):
							d[i,j] = ((np.sum(np.abs(df.iloc[i].values-self.get_entry(*str.split(self.closest_mill[self.V[i]],'_'))))) +
									  (np.sum(np.abs(self.get_entry(*str.split(self.closest_mill[self.V[i]],'_'))-df.iloc[j].values))))

		return d

	# Construct road distances matrix (returns both a distance and time matrix)
	def road_distance_matrix(self,df,mills = False):
		df = df.reset_index(drop=True).copy()
		d_dist = np.zeros([len(df),len(df)])
		d_time = np.zeros([len(df),len(df)])

		for i,row_i in df.iterrows():
			for j,row_j in df.iterrows():
				if i!=j:
					dist_time = road_distance.road_distance(df.iloc[i].values,df.iloc[j].values)
					d_dist[i,j] = dist_time['distance']
					d_time[i,j] = dist_time['time']
					if (self.M):
						if ((self.V[j] in self.H) and (self.V[i] in self.N)):
							seg_1 = road_distance.road_distance(df.iloc[i].values,self.get_entry(*str.split(self.closest_mill[self.V[i]],'_')))
							seg_2 = road_distance.road_distance(self.get_entry(*str.split(self.closest_mill[self.V[i]],'_')),df.iloc[j].values)
							d_dist[i,j] = seg_1['distance'] + seg_2['distance']
							d_time[i,j] = seg_1['time'] + seg_2['time']
		return (d_dist,d_time)



	def get_entry(self,a,v):
		if a == 'h':
			return self.H_p[int(v),:]
		if a == 'm':
			return self.M_p[int(v),:]
		if a == 'n':
			return self.N_p[int(v),:]

	# Calculate load violation
	def calc_load(self, route):
		return np.sum([self.quantities[n] for n in route])

	def total_load_violation(self, sol):
		return np.sum([np.max([self.calc_load(sol[h]['route']) - self.capacities[h],0]) for h in self.H])

	def path_to_format(self, path):
		new_path = []
		for n in path:
			if (n<len(N)):
				new_path.append("n_"+str(n))
			else:
				new_path.append("h_"+str(n-len(N)))
		return (new_path)

	def distance_path(self,  path):
		return np.sum([self.distance_dict[path[i]][path[i+1]] for i in range(len(path)-1)])

	def time_path(self,  path):
		return np.sum([self.time_dict[path[i]][path[i+1]] for i in range(len(path)-1)])

	def total_distance(self, sol):
		return np.sum([self.distance_path([h]+sol[h]['route']+[h]) for h in self.H])

	def set_sol_distance(self, sol):
		for h in self.H:
			sol[h]['route_dist'] = self.distance_path([h]+sol[h]['route']+[h])

	def LinearSolver(self, initialization = 0, time=False):
		H = self.H
		N = self.N
		H_p = self.H_p
		N_p = self.N_p
		quantities = self.quantities
		capacities = self.capacities
		type_dist = self.type_dist
		distance = self.distance

		def v_distance(v1, v2):
			return self.distance[self.positions[v1],self.positions[v2]]

		k = len(H)
		Y = ['y_'+ v for v in N]

		# Write model
		model = Model()

		#Create variables for the model
		vars = {}
		# First key is truck
		for truck in range(k):
			vars[truck] = {}
			# Going from house to farmer
			house = H[truck]
			vars[truck][house] = {}
			for farmer in N:
				vars[truck][house][farmer] = model.addVar(obj=v_distance(house, farmer), vtype=GRB.BINARY,
								  name='t_'+str(truck)+','+house+','+farmer)
			# Going from house to house
			vars[truck][house][house] = model.addVar(obj=0, vtype=GRB.BINARY,
								  name='t_'+str(truck)+','+house+','+house)
			# Going from farmer to farmer
			for f1 in N:
				vars[truck][f1] = {}
				for f2 in N:
					if (f1!=f2):
						vars[truck][f1][f2] = model.addVar(obj=v_distance(f1, f2), vtype=GRB.BINARY,
								  name='t_'+str(truck)+','+f1+','+f2)

			# Going from farmer to home
			for f in N:
				vars[truck][f][house] = model.addVar(obj=v_distance(f, house), vtype=GRB.BINARY,
								name='t_'+str(truck)+','+f+','+house)
		model.update()

		# We add restriction (1) so all farmers are visited
		for farmer in N:
			model.addConstr(quicksum(vars[truck][vertex][farmer] \
								 for truck in range(k) for vertex in (list(set(N)-set([farmer])) + [H[truck]]))==1)

		# Add restriction (2) so all nodes are visited once
		for truck in range(k):
			# This is for farmer nodes
			for farmer in N:
				s1 = quicksum(vars[truck][farmer][i] for i in vars[truck][farmer].keys())
				house = H[truck]
				s1 -= vars[truck][house][farmer]
				for f2 in N:
					if farmer!=f2:
						s1 -= vars[truck][f2][farmer]
				model.addConstr(s1 == 0)
			# This is for houses
			house = H[truck]
			s1 = quicksum(vars[truck][house][i] for i in vars[truck][house].keys())
			for farmer in N:
				s1 -= vars[truck][farmer][house]
			s1 -= vars[truck][house][house]
			model.addConstr(s1 == 0)


		# We add capacity constraints (3)
		for truck in range(k):
			s1 = 0
			for farmer in N:
				s1 += quantities[farmer]*quicksum(vars[truck][farmer][i] for i in vars[truck][farmer].keys())
			model.addConstr(s1 <= capacities['_'.join(['h',str(truck)])])

		# Trucks should leave their houses (4)
		for truck in range(k):
			house = H[truck]
			model.addConstr(quicksum(vars[truck][house][i] for i in vars[truck][house].keys()) == 1)

		# We finish with the cycle restrictions (5)
		y_vars = {v:model.addVar(obj=0, vtype=GRB.CONTINUOUS,
								  name='y_'+v) for v in (N)}
		model.update()

		for i in y_vars.keys():
			for j in y_vars.keys():
				if i!=j:
					s1 = y_vars[i] - y_vars[j]
					for truck in range(k):
						if j in vars[truck][i].keys():
							s1 += len(y_vars) * vars[truck][i][j]
					model.addConstr(s1 <= len(y_vars) - 1)

		model.update()

		model.reset()
		model._vars = model.getVars()

		if initialization:
			# We create the list of variables that are positive
			positive = []
			for h in initialization.keys():
				path = [h]+initialization[h]['route']+[h]
				t = '_'.join(['t',h.split('_')[1]])
				for i in range(len(path)-1):
					positive.append(','.join([t,path[i],path[(i+1)]]))
			y_init = {}
			i = 0
			for h in initialization.keys():
				path = initialization[h]['route']
				for s in path:
					if s[0] == 'n':
						y_init['y_'+s] = i
						i = i+1


			# Initialize variables
			variables = model.getVars()
			for j,v in enumerate(variables):
				name = v.VarName
				if name[0]=='t':
					if name in positive:
						model.getVars()[j].start = 1.0
					else:
						model.getVars()[j].start = 0.0
				elif name[0]=='y':
					model.getVars()[j].start = y_init[name]

		if time:
			model.setParam('TimeLimit', time)
		model.optimize()
		self.model = model
		solution = self.solution_from_solver(model)

		# Edit solution
		for h in H:
			l = solution[h]['route']
			solution[h]['load'] = np.sum([quantities[n] for n in l])
			solution[h]['max_load'] = capacities[h]

		self.set_sol_distance(solution)

		return solution

	# Takes a linear model, returns the solution of that model
	def solution_from_solver(self, model):
		# Reconstruct the solution
		solution_dict = {}
		for v in model.getVars():
			name = v.VarName
			value = v.x
			if value == 1.0 and name[0] == 't':
				parts = name.split(',')
				truck = parts[0]
				node_0 = parts[1]
				node_1 = parts[2]
				if not truck in solution_dict.keys():
					solution_dict[truck] = {}
				solution_dict[truck][node_0] = node_1
		sol_from_model = {}
		for t in solution_dict.keys():
			h = '_'.join(['h',t.split('_')[1]])
			sequence = []
			k = h
			while True:
				k = solution_dict[t][k]
				if k==h:
					break
				sequence.append(k)
			sol_from_model[h] = {'route':sequence}
		return sol_from_model

	def draw_solution(self, solution):
		# Draw the solution of the heuristic
		H = self.H
		N = self.N
		H_p = self.H_p
		N_p = self.N_p
		quantities = self.quantities
		capacities = self.capacities
		type_dist = self.type_dist
		M_p = self.M_p
		M = self.M


		plt.scatter(N_p[:,0], N_p[:,1], c='r')
		plt.scatter(H_p[:,0], H_p[:,1], c='r', marker = 's', s = 100)
		if self.M:
			plt.scatter(M_p[:,0], M_p[:,1], c='r', marker = 'o', s = 100)
		for truck in solution.keys():
			sol = solution[truck]
			route = sol['route']
			if len(route) > 0:
				for i in range(len(route) - 1):
					coord1 = self.get_entry(*str.split(route[i],'_'))
					coord2 = self.get_entry(*str.split(route[i+1],'_'))
					plt.plot([coord1[0], coord2[0]], [coord1[1], coord2[1]], color='r', linestyle='-', linewidth=1)
				coord1 = self.get_entry(*str.split(truck,'_'))
				coord2 = self.get_entry(*str.split(route[0],'_'))
				plt.plot([coord1[0], coord2[0]], [coord1[1], coord2[1]], color='r', linestyle='-', linewidth=1)
				coord1 = self.get_entry(*str.split(route[len(route) - 1],'_'))
				coord2 = self.get_entry(*str.split(truck,'_'))
				if self.M:
					coord3 = self.get_entry(*str.split(self.closest_mill[route[len(route) - 1]],'_'))
					plt.plot([coord1[0], coord3[0]], [coord1[1], coord3[1]], color='r', linestyle='-', linewidth=1)
					plt.plot([coord3[0], coord2[0]], [coord3[1], coord2[1]], color='r', linestyle='-', linewidth=1)
				else:
					plt.plot([coord1[0], coord2[0]], [coord1[1], coord2[1]], color='r', linestyle='-', linewidth=1)
		plt.show()

	def draw_problem(self):
		# Draw the solution of the heuristic
		H = self.H
		N = self.N
		H_p = self.H_p
		N_p = self.N_p
		quantities = self.quantities
		capacities = self.capacities
		type_dist = self.type_dist
		M_p = self.M_p
		M = self.M


		plt.scatter(N_p[:,0], N_p[:,1], c='r')
		plt.scatter(H_p[:,0], H_p[:,1], c='r', marker = 's', s = 100)
		if self.M:
			plt.scatter(M_p[:,0], M_p[:,1], c='r', marker = 'o', s = 100)
		plt.show()

	def draw_path(self, path):
		# Draw the solution of the heuristic
		H = self.H
		N = self.N
		H_p = self.H_p
		N_p = self.N_p
		quantities = self.quantities
		capacities = self.capacities
		type_dist = self.type_dist
		M_p = self.M_p
		M = self.M


		plt.scatter(N_p[:,0], N_p[:,1], c='r')
		plt.scatter(H_p[:,0], H_p[:,1], c='r', marker = 's', s = 100)
		if self.M:
			plt.scatter(M_p[:,0], M_p[:,1], c='r', marker = 'o', s = 100)
		route = path
		if len(route) > 0:
			for i in range(len(route) - 1):
				coord1 = self.get_entry(*str.split(route[i],'_'))
				coord2 = self.get_entry(*str.split(route[i+1],'_'))
				plt.plot([coord1[0], coord2[0]], [coord1[1], coord2[1]], color='r', linestyle='-', linewidth=1)
		plt.show()

	def random_start(self):

		H = self.H
		N = self.N
		H_p = self.H_p
		N_p = self.N_p
		quantities = self.quantities
		capacities = self.capacities
		distance = self.distance
		distance_dict = self.distance_dict

		assignment_dict = {h:{'route':[],'cap':0} for h in H}
		for n in N:
			dist_to_trucks = [np.random.rand() for h in H]
			ranking = np.argsort(dist_to_trucks)
			q = quantities[n]
			assign = False
			for r in ranking:
				truck = H[r]
				truck_capacity = capacities[truck]
				current_load = assignment_dict[truck]['cap']
				if current_load + q <= truck_capacity:
					l = assignment_dict[truck]['route']
					l.append(n)
					assignment_dict[truck]['cap'] = current_load + q
					assignment_dict[truck]['route'] = l
					assign = True
					break
			if assign == False:
				print("Problem infeasible")

		# Edit solution
		for h in H:
			l = assignment_dict[h]['route']
			assignment_dict[h]['load'] = np.sum([quantities[n] for n in l])
			assignment_dict[h]['max_load'] = capacities[h]

		self.set_sol_distance(assignment_dict)

		return assignment_dict

	def start_heuristic(self):
		H = self.H
		N = self.N
		H_p = self.H_p
		N_p = self.N_p
		quantities = self.quantities
		capacities = self.capacities
		distance = self.distance
		distance_dict = self.distance_dict

		assignment_dict = {h:{'nodes':[],'cap':0} for h in H}
		for n in N:
			dist_to_trucks = [distance_dict[n][h] for h in H]
			ranking = np.argsort(dist_to_trucks)
			q = quantities[n]
			assign = False
			for r in ranking:
				truck = H[r]
				truck_capacity = capacities[truck]
				current_load = assignment_dict[truck]['cap']
				if current_load + q <= truck_capacity:
					l = assignment_dict[truck]['nodes']
					l.append(n)
					assignment_dict[truck]['cap'] = current_load + q
					assignment_dict[truck]['nodes'] = l
					assign = True
					break
			if assign == False:
				print("Problem infeasible")
		assignment_nodes = {n:h  for h in assignment_dict.keys() for n in assignment_dict[h]['nodes']}

		# Create the routes of the heuristic
		iterations = 100
		solution_tsp = {}
		for h in H:
			if len(assignment_dict[h]['nodes'])>0:
				nodes = set(assignment_dict[h]['nodes'] + [h])
				best_dist = float('inf')
				best_route = list(nodes)
				for i in range(iterations):
					nn_solver = NearestNeighbourSolver(nodes, distance_dict)
					new_route = nn_solver.run()
					new_dist = nn_solver.get_total_distance()
					if new_dist < best_dist:
						best_dist = new_dist
						best_route = new_route
					ni_solver = NearestInsertionSolver(nodes, distance_dict)
					new_route = ni_solver.run()
					new_dist = ni_solver.get_total_distance()
					if new_dist < best_dist:
						best_dist = new_dist
						best_route = new_route
				solution_tsp[h] = {'route':best_route}
			else:
				solution_tsp[h] = {'route':[]}

		# Edit solution
		for h in H:
			l = solution_tsp[h]['route']
			if len(solution_tsp[h]['route'])>0:
				l.pop()
				ind = l.index(h)
				l = l[ind+1:len(l)]+l[0:ind]
				solution_tsp[h]['route'] = l
			solution_tsp[h]['load'] = np.sum([quantities[n] for n in l])
			solution_tsp[h]['max_load'] = capacities[h]

		self.set_sol_distance(solution_tsp)

		return solution_tsp

	def heuristic(self, solution):

		H = self.H
		N = self.N
		H_p = self.H_p
		N_p = self.N_p
		quantities = self.quantities
		capacities = self.capacities
		type_dist = self.type_dist
		distance = self.distance
		distance_dict = self.distance_dict

		solution_tsp = copy.deepcopy(solution)

		attributes = {}
		for truck in H:
			for node in N:
				attributes[(truck,node)] = {'tabu':0,'rho':0}

		def rho(sol):
			return np.sum([ attributes[(h,n)]['rho'] for h in H for n in sol[h]['route']])

		alpha = 0.015
		lamb = 1
		theta = int(7.5*np.log10(len(N)))
		delta = 0.5
		# Calculate current cost
		dist = self.total_distance(solution_tsp)
		load = self.total_load_violation(solution_tsp)
		f = dist + alpha*load
		rho_sum = rho(solution_tsp)

		best_dist = dist
		best_sol = copy.deepcopy(solution_tsp)
		for it in range(1000):
			partial_data = {
				'cost' : float('inf'),
				'truck' : -1,
				'new_truck' : -1,
				'index' : -1,
				'new_index' : -1,
				'node' : -1
			}

			# Find two attributes to change
			for truck in H:
				for node in solution_tsp[truck]['route']:
					complete_route = [truck]+solution_tsp[truck]['route']+[truck]
					ind_node = complete_route.index(node)

					dec_dist = (distance_dict[complete_route[ind_node-1]][complete_route[ind_node+1]]
								- distance_dict[complete_route[ind_node-1]][complete_route[ind_node]]
								- distance_dict[complete_route[ind_node+1]][complete_route[ind_node]])
					overload = np.max([solution_tsp[truck]['load'] - solution_tsp[truck]['max_load'],0])
					dec_overload = - np.min([overload, quantities[node]])

					for new_truck in H:
						if truck != new_truck:
							overload = np.max([solution_tsp[new_truck]['load'] - solution_tsp[new_truck]['max_load'],0])
							new_overload = np.max([solution_tsp[new_truck]['load'] + quantities[node] - solution_tsp[new_truck]['max_load'],0])
							inc_overload_truck = new_overload - overload
							change_overload = inc_overload_truck + dec_overload
							new_complete_route = [new_truck]+solution_tsp[new_truck]['route']+[new_truck]
							for i in range(len(new_complete_route)-1):
								inc = (distance_dict[new_complete_route[i]][node]
									   + distance_dict[new_complete_route[i+1]][node]
									   -distance_dict[new_complete_route[i]][new_complete_route[i+1]])
								change_dist = inc + dec_dist

								# Calculate the new cost
								new_distance = (dist + change_dist)
								new_load = (load + change_overload)
								new_f = new_distance + alpha*new_load
								new_rho = rho_sum - attributes[(truck,node)]['rho'] + attributes[(new_truck,node)]['rho']
								new_penalty = lamb*new_distance*np.sqrt(len(N)*len(H))*new_rho
								cost = new_f+new_penalty

								if (cost < partial_data['cost']) and (attributes[(new_truck,node)]['tabu']==0 or new_f<f) :
									partial_data = {
										'cost' : cost,
										'truck' : truck,
										'new_truck' : new_truck,
										'index' : ind_node,
										'new_index' : i,
										'node' : node,
										'dist_gain' : change_dist,
										'load_gain' : change_overload
									}

			# Do the update of the solution
			old_path = solution_tsp[partial_data['truck']]['route']
			old_path.pop(partial_data['index'] - 1)
			solution_tsp[partial_data['truck']]['route'] = old_path
			new_path = solution_tsp[partial_data['new_truck']]['route']
			new_path.insert(partial_data['new_index'],partial_data['node'])
			solution_tsp[partial_data['new_truck']]['route'] = new_path

			# Update loads
			solution_tsp[partial_data['new_truck']]['load'] += quantities[partial_data['node']]
			solution_tsp[partial_data['truck']]['load'] -= quantities[partial_data['node']]

			# Update attributes
			attributes[(partial_data['new_truck'],partial_data['node'])]['rho'] += 1
			for a in attributes.keys():
				attributes[a]['tabu'] = np.max([attributes[a]['tabu']-1,0])
			attributes[(partial_data['truck'],partial_data['node'])]['tabu'] = theta
			# Also update global vars
			dist = self.total_distance(solution_tsp)
			load = self.total_load_violation(solution_tsp)
			rho_sum = rho(solution_tsp)
			f = dist + alpha*load

			if load>0:
				alpha = alpha * (1+delta)
			else:
				alpha = alpha / (1+delta)
				if dist<best_dist:
					best_dist = dist
					best_sol = copy.deepcopy(solution_tsp)
					print(it, best_dist)
		self.set_sol_distance(best_sol)
		return best_sol

	def christofides_bound(self, iterations, z_ub, epsilon):
		H = self.H
		N = self.N
		# We need to discretize the quantities
		# Run the first procedure
		distance_dictionary = copy.deepcopy(self.distance_dict)
		int_quantities = {n:round(self.quantities[n]*10) for n in self.quantities.keys()}
		int_capacities = {h:round(self.capacities[h]*10) for h in self.capacities.keys()}

		max_val,u_opt,v_opt,lamb_opt = lower_bound.optimize_lower_bound(iterations, z_ub, epsilon, H,int_capacities,N,int_quantities,distance_dictionary)
		return (max_val)

	def feasibility(self):
		model = Model()
		capacities = self.capacities
		quantities = self.quantities

		variables = {}
		for t, cap in capacities.iteritems():
			variables[t] = {}
			for f, q in quantities.iteritems():
				variables[t][f] = model.addVar(obj=0, vtype=GRB.BINARY,
											  name=t+f)
		model.update()
		for t, cap in capacities.iteritems():
			model.addConstr(quicksum([variables[t][f] * quantities[f]  for f in quantities.keys()])<=capacities[t])

		for f, q in quantities.iteritems():
			model.addConstr(quicksum([variables[t][f]  for t in capacities.keys()])==1)
		model.setParam('OutputFlag', 0)
		model.update()
		model.optimize()
		model.feasibility()
		if model.Status == 3:
			return(False)
		else:
			return(True)

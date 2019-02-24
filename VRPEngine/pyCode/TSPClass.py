import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import math
import random
from gurobipy import *
import geopy.distance
import pandas as pd
import copy
from scipy.spatial import distance as distance_lib
import road_distance
import sys
sys.path.insert(0,'/Users/sergiocamelo/Dropbox/Sergio-Joann/Code/VRPEngine/pyCode/pyconcorde')
from concorde.tsp import TSPSolver
from concorde.tests.data_utils import get_dataset_path

# This class takes E,N (set of extremes and nodes), as well as their positions. A distance matrix that has first
# the extremes and then the nodes is passed
class TSP:
	def __init__(self, E, N, E_p, N_p, type_dist, distance_matrix = None):
		self.E = E
		self.N = N
		self.V = E + N
		self.E_p = E_p
		self.N_p = N_p
		self.type_dist = type_dist

		if distance_matrix:
			self.distance = distance_matrix
		else:
			if type_dist == 'road':
				aux = self.road_distance_matrix(pd.DataFrame(data = np.concatenate([E_p, N_p]), index = np.concatenate([self.E,self.N]), columns = ['lat','lon']))
				self.road_distance = aux[0]
				self.road_time = aux[1]
				self.distance = self.road_distance
				self.time = self.road_time

			if type_dist == 'euclid':
				self.distance = self.euclid_distance_matrix(pd.DataFrame(data = np.concatenate([E_p, N_p]), index = np.concatenate([self.E,self.N]), columns = ['lat','lon']))



	def solve(self):

		if (len(self.N) <= 1):
			cycle = self.E+self.N
			return {
			"solution":cycle,
			"objective":np.sum([self.distance[i][i+1] for i in range(len(cycle)-1)])
			}


		precision = 1
		# Let's construct the file
		tsp_template = """TYPE : TSP
DIMENSION : {}
EDGE_WEIGHT_TYPE : EXPLICIT
EDGE_WEIGHT_FORMAT : LOWER_DIAG_ROW 
EDGE_WEIGHT_SECTION :
{}
EOF"""
		len_v = len(self.V)
		distances_data = ''
		for i in range(0,len_v):
			for j in range(0,i+1):
				if i == j: 
					distances_data += '0\n'
				else: 
					distances_data += '{} '.format(int(self.distance[i][j] * precision))

		# Create a dummy variable
		for j in range(0,len_v+1):
			if j == 0 or j == 1 or j == len_v:
				distances_data += '{} '.format(0)
			else:
				distances_data += '{} '.format(int(np.max(self.distance * precision) * 10e3))


		tsp_data = tsp_template.format(len_v+1, distances_data)
		print(tsp_data)
		with open("/Users/sergiocamelo/Dropbox/Sergio-Joann/tsp.TSPLIB", "w") as text_file:
			text_file.write(tsp_data)
		solver = TSPSolver.from_tspfile("/Users/sergiocamelo/Dropbox/Sergio-Joann/tsp.TSPLIB")

		solution = solver.solve()
		tour = solution.tour
		objective = solution.optimal_value/(precision+0.0)
		print(tour)

		cycle = []
		dummy_position = list(tour).index(len_v)
		for i in range(dummy_position+1,len_v+1):
			cycle.append(self.V[tour[i]])
		for i in range(0,dummy_position):
			cycle.append(self.V[tour[i]])
		return {
			"solution":cycle,
			"objective":objective
		}




	# Construct geographical distances matrix
	def geo_distance_matrix(self,df):
		df = df.reset_index(drop=True).copy()
		d = np.zeros([len(df),len(df)])
		for i,row_i in df.iterrows():
			for j,row_j in df.iterrows():
				if i!=j:
					d[i,j] = geopy.distance.vincenty(df.iloc[i].values,df.iloc[j].values).m
		return d

	# Construct euclid distances matrix
	def euclid_distance_matrix(self,df):
		df = df.reset_index(drop=True).copy()
		d = np.zeros([len(df),len(df)])
		for i,row_i in df.iterrows():
			for j,row_j in df.iterrows():
				if i!=j:
					d[i,j] = np.sqrt(np.sum((df.iloc[i].values-df.iloc[j].values)**2))

		return d

	# Construct manhattan distances matrix
	def manhattan_distance_matrix(self,df):
		df = df.reset_index(drop=True).copy()
		d = np.zeros([len(df),len(df)])
		for i,row_i in df.iterrows():
			for j,row_j in df.iterrows():
				if i!=j:
					d[i,j] = np.sum(np.abs(df.iloc[i].values-df.iloc[j].values))

		return d

	# Construct road distances matrix (returns both a distance and time matrix)
	def road_distance_matrix(self,df):
		df = df.reset_index(drop=True).copy()
		d_dist = np.zeros([len(df),len(df)])
		d_time = np.zeros([len(df),len(df)])

		for i,row_i in df.iterrows():
			for j,row_j in df.iterrows():
				if i!=j:
					dist_time = road_distance.road_distance(df.iloc[i].values,df.iloc[j].values)
					d_dist[i,j] = dist_time['distance']
					d_time[i,j] = dist_time['time']
		return (d_dist,d_time)

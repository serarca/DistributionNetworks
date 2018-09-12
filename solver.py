# This script gathers the solvers for our optimization problem
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist
import math
import random
from gurobipy import *
import geopy.distance

# Solves capacitated VRP, receives arrays of nodes, arrays of positions, dictionaries of quantities and capacities,
# a distance matrix in the shape of a Pandas array indexed by the nodes, and a boolean 
def LinearSolver(H, M, N, H_p, M_p, N_p, quantities, capacities, distance, euclidean):
    values = {'H':H_p,'M':M_p,'N':N_p}
    # Auxiliary function
    def get_entry(a,v):
        t=v
        if a == 'h':
            return values['H'][int(v),:]
        if a == 'm':
            return values['M'][int(v),:]
        if a == 'n':
            return values['N'][int(v),:]
    # These functions are auxiliary in the constructin of the distance matrix
    def v_distance(v1, v2):
        if euclidean:
            return geopy.distance.vincenty(get_entry(*str.split(v1,'_')),get_entry(*str.split(v2,'_'))).m
        else:
            return distance.loc[v1,v2]

    k = len(H)
    Y = ['y_'+ v for v in M + N]
    
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

        # Going from farmer to mill
        for farmer in N:
            for mill in M:
                vars[truck][farmer][mill] = model.addVar(obj=v_distance(farmer, mill), vtype=GRB.BINARY,
                            name='t_'+str(truck)+','+farmer+','+mill)

        # Going from mill to home
        for mill in M:
            vars[truck][mill] = {}
            vars[truck][mill][house] = model.addVar(obj=v_distance(mill, house), vtype=GRB.BINARY,
                            name='t_'+str(truck)+','+mill+','+house)
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
        for m in M:
            s1 -= vars[truck][m][house]
        s1 -= vars[truck][house][house]
        model.addConstr(s1 == 0)
        # This is for mills
        for mill in M:
            s1 = quicksum(vars[truck][mill][i] for i in vars[truck][mill].keys())
            for farmer in N:
                s1 -= vars[truck][farmer][mill]
            model.addConstr(s1 == 0)

    # We add capacity constraints (3)
    for truck in range(k):
        s1 = 0
        for farmer in N:
            s1 += quantities[farmer]*quicksum(vars[truck][farmer][i] for i in vars[truck][farmer].keys())
        model.addConstr(s1 <= capacities[truck])

    # Trucks should leave their houses (4)
    for truck in range(k):
        house = H[truck]
        model.addConstr(quicksum(vars[truck][house][i] for i in vars[truck][house].keys()) == 1)

    # We finish with the cycle restrictions (5)
    y_vars = {v:model.addVar(obj=0, vtype=GRB.CONTINUOUS,
                              name='y_'+v) for v in (M + N)} 
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
    #model.optimize(mycallback)

    return model

# Draw points
# Mills are in yellow circles
# Houses of trucks are squares
# Dots are farmers
def draw_problem(H_p, M_p, N_p, quantities_v):
    colors = ['r', 'b', 'c', 'm']
    plt.scatter(N_p[:,0], N_p[:,1], c='k', s = quantities_v)
    plt.scatter(M_p[:,0], M_p[:,1], c='y', marker = 'o', s = 200)
    plt.scatter(H_p[:,0], H_p[:,1], c= colors[0:len(H_p)], marker = 's', s = 100)
    plt.show()
    
def draw_solution(H_p, M_p, N_p, quantities_v, model):
    colors = ['r', 'b', 'c', 'm','r', 'b', 'c', 'm','r', 'b', 'c', 'm','r', 'b', 'c', 'm']
    values = {'H':H_p,'M':M_p,'N':N_p}
    def get_entry(a,v):
        t=v
        if a == 'h':
            return values['H'][int(v),:]
        if a == 'm':
            return values['M'][int(v),:]
        if a == 'n':
            return values['N'][int(v),:]
    plt.scatter(N_p[:,0], N_p[:,1], c='k', s = quantities_v)
    plt.scatter(M_p[:,0], M_p[:,1], c='y', marker = 'o', s = 200)
    plt.scatter(H_p[:,0], H_p[:,1], c= colors[0:len(H_p)], marker = 's', s = 100)
    variables = model.getVars()
    for v in variables:
        if v.varName[0] == 't' and v.x == 1.0:
            names = v.varName.split(',')
            p1 = get_entry(*str.split(names[1],'_'))
            x1 = p1[0]
            y1 = p1[1]
            p2 = get_entry(*str.split(names[2],'_'))
            x2 = p2[0]
            y2 = p2[1]
            c = colors[int(str.split(names[0],'_')[1])]
            plt.plot([x1, x2], [y1, y2], color=c, linestyle='-', linewidth=1)
    plt.show()
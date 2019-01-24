import random
import sys
sys.path.insert(0,'/Users/sergiocamelo/Dropbox/Sergio-Joann/Code/VRPEngine/pyCode/')
from TSPClass import TSP


n = 10
N = ['n_'+str(i) for i in range(n)]
E = ['e_0','e_1']

E_p = [[0,0],[1,1]]
N_p = [[random.random() for i in range(2)] for j in range(n)]


tsp = TSP(E, N, E_p, N_p, 'euclid')
tsp.solve()
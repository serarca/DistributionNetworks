import networkx as nx

# Create Network class
class NetworkClass:
    def __init__(self, network):
        self.network = network
        
    def network_distance(self,node1,node2,weight):        
        return nx.shortest_path_length(self.network, node1, node2, weight=weight)
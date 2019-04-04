# This class holds a VRP problem
class PrimalSolution:
	def __init__(self, routes, z_lb):
		self.routes = routes
		self.z_lb = z_lb
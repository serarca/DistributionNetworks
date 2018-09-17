import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.string cimport string
from libcpp.list cimport list
from libcpp.set cimport set



cdef extern from "lower_bounds.h":
   struct QPaths:
      vector[vector[double]] f
      vector[vector[double]] phi
      vector[vector[int]] p
      vector[vector[vector[int]]] q_route
      vector[vector[vector[int]]] q_route_2
   struct QRoutes:
      vector[vector[double]] psi
      vector[vector[vector[int]]] psi_route
   struct LowerBound:
      double z_lb
      vector[double] theta
      vector[double] rho
      vector[double] u
   struct DualSolution:
      double z_lb
      vector[double] lamb
      vector[double] u
      vector[double] v
      vector[list[SimpleRoute]] routes
      double gamma_guarantee

   QPaths construct_q_paths_(
      int h,
      int truck_capacity,
      vector[int] N,
      vector[vector[double]] distance,
      vector[int] values,
      map[int,int] values_pos,
      vector[int] quantities,
      string direction
      )
   QRoutes construct_q_routes_(
      int h,
      int truck_capacity,
      vector[int] N,
      vector[vector[double]] distance,
      vector[int] values,
      map[int,int] values_pos,
      vector[int] quantities
      )

   LowerBound lower_bound_(
      vector[int] H,
      vector[int] capacities,
      vector[int] N,
      vector[int] quantities,
      vector[vector[double]] distance_dict,
      vector[double] mu,
      vector[double] lamb,
      vector[int] n_trucks
      )
   DualSolution lower_bound_optimizer_M1(
      int iterations,
      double z_ub,
      double epsilon,
      vector[int] H,
      vector[int] capacities,
      vector[int] N,
      vector[int] quantities,
      vector[vector[double]] geo_distance,
      vector[int] n_trucks);


cdef extern from "baldacci.h":
   struct Path:
      list[int] path
      set[int] nodes
      double cost
      double lower_bound
      int load
      int end

   struct SimpleRoute:
      list[int] path
      int index_l
      int index_r
      double cost
      int load
      int median
      double geo_cost
      int truck
      double l_lb
      double r_lb


   vector[DualSolution] construct_lower_bound(
      int iterations_grad_m1,
      int iterations_grad_m2,
      int iterations_m2,
      double z_ub,
      int Delta,
      int Delta_zero,
      int Delta_final,
      double gamma,
      double gamma_zero,
      double gamma_final,
      double epsilon,
      vector[int] H,
      vector[int] capacities,
      vector[int] N,
      vector[int] quantities,
      vector[vector[double]] geo_distance,
      vector[int] n_trucks
   )

cpdef construct_q_paths(h_,truck_capacity_,N_,distance_,values_,values_pos_,quantities_,direction_):
    cdef:
        int h = h_
        int truck_capacity = truck_capacity_
        vector[int] N = N_
        vector[vector[double]] distance = distance_
        vector[int] values = values_
        map[int,int] values_pos = values_pos_
        vector[int] quantities = quantities_
        string direction = direction_

    cdef QPaths qpaths = construct_q_paths_(h,truck_capacity,N,distance,values,values_pos,quantities,direction)
    return qpaths

cpdef construct_q_routes(h_,truck_capacity_,N_,distance_,values_,values_pos_,quantities_):
   cdef:
      int h = h_
      int truck_capacity = truck_capacity_
      vector[int] N = N_
      vector[vector[double]] distance = distance_
      vector[int] values = values_
      map[int,int] values_pos = values_pos_
      vector[int] quantities = quantities_

   cdef QRoutes qroutes = construct_q_routes_(h,truck_capacity,N,distance,values,values_pos,quantities)
   return qroutes

cpdef lower_bound(H_, capacities_, N_, quantities_, distance_, mu_, lamb_, n_trucks_):
   cdef:
      vector[int] H = H_,
      vector[int] capacities = capacities_,
      vector[int] N = N_,
      vector[int] quantities = quantities_,
      vector[vector[double]] distance = distance_,
      vector[double] mu = mu_,
      vector[double] lamb = lamb_
      vector[int] n_trucks = n_trucks_
   cdef LowerBound lb = lower_bound_(H, capacities, N, quantities, distance, mu, lamb, n_trucks)
   return lb


cpdef lower_bound_optimizer_M1_(iterations_, z_ub_, epsilon_, H_, capacities_, N_, quantities_, geo_distance_, n_trucks_):
   cdef:
      int iterations = iterations_
      double z_ub = z_ub_
      double epsilon = epsilon_
      vector[int] H = H_
      vector[int] capacities = capacities_
      vector[int] N = N_
      vector[int] quantities = quantities_
      vector[vector[double]] geo_distance = geo_distance_
      vector[int] n_trucks = n_trucks_
   cdef DualSolution dual_solution = lower_bound_optimizer_M1(iterations, z_ub, epsilon, H, capacities, N, quantities, geo_distance, n_trucks)
   return dual_solution

cpdef construct_lower_bound_(iterations_grad_m1_,iterations_grad_m2_,iterations_m2_,z_ub_,Delta_,Delta_zero_,Delta_final_,gamma_,gamma_zero_,gamma_final_,epsilon_,H_,capacities_,N_,quantities_,geo_distance_,n_trucks_):
   cdef:
      int iterations_grad_m1 = iterations_grad_m1_
      int iterations_grad_m2 = iterations_grad_m2_
      int iterations_m2 = iterations_m2_
      double z_ub = z_ub_
      int Delta = Delta_
      int Delta_zero = Delta_zero_
      int Delta_final = Delta_final_
      double gamma = gamma_
      double gamma_zero = gamma_zero_
      double gamma_final = gamma_final_
      double epsilon = epsilon_
      vector[int] H = H_
      vector[int] capacities = capacities_
      vector[int] N = N_
      vector[int] quantities = quantities_
      vector[vector[double]] geo_distance = geo_distance_
      vector[int] n_trucks = n_trucks_
   cdef vector[DualSolution] dual_solution = construct_lower_bound(iterations_grad_m1,iterations_grad_m2,iterations_m2,z_ub,Delta,Delta_zero,Delta_final,gamma,gamma_zero,gamma_final,epsilon,H,capacities,N,quantities,geo_distance,n_trucks)
   return dual_solution

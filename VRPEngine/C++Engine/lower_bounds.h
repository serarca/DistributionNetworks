#pragma once


#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <list>
#include "VRPClass.cpp"



using namespace std;



// The struct of results for q-routes
struct QRoutes {
   vector<vector<double>> psi;
   vector<vector<vector<int>>> psi_route;
};



// The struct of results
struct QPaths {
   vector<vector<double>> f;
   vector<vector<double>> phi;
   vector<vector<int>> p;
   vector<vector<vector<int>>> q_route;
   vector<vector<vector<int>>> q_route_2;
};

// The struct of lower bounds
struct LowerBound {
   double z_lb;
   vector<double> theta;
   vector<double> rho;
   vector<double> u;
};

//The struct of paths
struct SimpleRoute{
   list<int> path;
   int index_l;
   int index_r;
   double cost;
   int load;
   int median;
   double geo_cost;
   int truck;
   double l_lb;
   double r_lb;
};


class DualSolution {
public:
   double z_lb;
   vector<double> lamb;
   vector<double> u;
   vector<double> v;
   vector<list<SimpleRoute>> routes;
   double gamma_guarantee;
   vector<vector<double>> reduced_distances;

   void calc_reduced_distances(vector<vector<double>> geo_distance);
   void initialize_routes(int len_H);

};

struct TerminatingCondition {
   bool terminated;
   double gamma_guarantee;
   int new_routes;
};



struct LB_GENPATH {
   vector<vector<double>> F;
   vector<vector<double>> G;
   vector<vector<int>> X;
   vector<vector<vector<int>>> min_q_path;
   vector<vector<vector<int>>> min_q_path_2;
};

vector<vector<double>> reduced_cost_matrix(
   vector<vector<double>> geo_distance,
   vector<double> lamb,
   vector<double> mu
);

PossibleValues possible_values(vector<int>& quantities, int truck_capacity);

LowerBound lower_bound_(
   VRP &vrp,
   vector<vector<double>> &distance_dict,
   vector<double> mu,
   vector<double> lamb
);

QPaths construct_q_paths_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> &distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities,
   string direction
);

QRoutes construct_q_routes_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> &distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities
);

DualSolution lower_bound_optimizer_M1(
   int iterations,
   double z_ub,
   double epsilon,
   VRP &vrp
);

LB_GENPATH path_lower_bound(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> &distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities,
   string direction
);

vector<vector<double>> penalized_matrix(
   vector<vector<double>> distance_dict,
   vector<double> penalties,
   int len_N,
   int len_H,
   double penalty_factor
);

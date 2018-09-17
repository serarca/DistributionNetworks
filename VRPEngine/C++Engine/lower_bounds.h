#pragma once


#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <list>


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


struct DualSolution {
   double z_lb;
   vector<double> lamb;
   vector<double> u;
   vector<double> v;
   vector<list<SimpleRoute>> routes;
   double gamma_guarantee;
};

// A struct with posible values and their inverse maps
struct PossibleValues {
   vector<int> values;
   map<int,int> values_pos;
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
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> distance_dict,
   vector<double> mu,
   vector<double> lamb,
   vector<int> n_trucks
);

QPaths construct_q_paths_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities,
   string direction
);

QRoutes construct_q_routes_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities
);

DualSolution lower_bound_optimizer_M1(
   int iterations,
   double z_ub,
   double epsilon,
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> geo_distance,
   vector<int> n_trucks
);

LB_GENPATH path_lower_bound(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities,
   string direction
);

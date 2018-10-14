#pragma once

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <list>
#include <set>
#include "lower_bounds.h"
#include <bitset>

using namespace std;

typedef std::bitset<100> bit_set;


//The struct of paths
struct Path{
   list<int> path;
   bit_set nodes;
   double cost;
   double lower_bound;
   int load;
   int end;
};

//The struct of paths
struct Route{
   list<Path>::iterator path_l;
   list<Path>::iterator path_r;
   int index_l;
   int index_r;
   bit_set nodes;
   double cost;
   int load;
   int median;
   double l_lb;
   double r_lb;
};




vector<list<Path>> GENPATH(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> &distance_dict,
   string direction,
   bool &terminated,
   double &gamma_guarantee
);

list<SimpleRoute> GENROUTE(
   double z_ub,
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> &distance_dict,
   vector<vector<double>> &geo_distance,
   bool &terminated,
   double &gamma_guarantee
);

DualSolution optimize_lower_bound_M2(
   int sub_iterations,
   double z_ub,
   int Delta,
   int Delta_zero,
   double gamma,
   double gamma_zero,
   double epsilon,
   VRP &vrp,
   vector<double> mu,
   vector<double> lamb,
   vector<double> u,
   vector<list<SimpleRoute>> &initial_routes,
   double &initial_gamma_guarantee
);

vector<DualSolution> construct_lower_bound(
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
   VRP &vrp
);

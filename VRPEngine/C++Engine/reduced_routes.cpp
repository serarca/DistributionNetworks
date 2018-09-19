#pragma once


#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <list>
#include "lower_bounds.h"
#include "baldacci.h"

#include<cmath>
#include "prettyprint.hpp"
#include <bitset>


vector<list<SimpleRoute>> get_reduced_routes(
   double z_ub,
   int Delta,
   double gamma,
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> geo_distance,
   vector<double> mu,
   vector<double> u
){

	  //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();

   // Define the vector of routes for each truck
   vector<list<SimpleRoute>> Routes(len_H);

   // Calculate reduced costs
   vector<vector<double>> distance_dict = reduced_cost_matrix(geo_distance, u, mu);

   // Define infinity
   double inf = numeric_limits<double>::infinity();

   // We start by generating routes for all of the trucks
   bool terminated_initial = true;
   for (auto h:H){
      cout<<"Generating Routes for Truck: "<<h<<endl;
      double truck_guarantee;
      Routes[h - len_N] = GENROUTE(z_ub, Delta, gamma, h, capacities[h - len_N], N, quantities, distance_dict, geo_distance, terminated_initial, truck_guarantee);
   }
   return Routes;
}





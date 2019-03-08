#pragma once

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <numeric>
#include "lower_bounds.h"
#include <algorithm>
#include <numeric>
#include <math.h>
#include "prettyprint.hpp"
#include "VRPClass.cpp"
#include "lower_bounds.cpp"



// This function returns the q-paths that the q-routes function uses
// to calulate lower bounds
QPaths construct_n_q_paths_(
   VRP& vrp,
   int h,
   string direction,
   int limit
){

   //Calculate lengths
   int len_values = vrp.possible_v[h-vrp.len_N()].values.size();
   int len_N = vrp.len_N();

   //Construct infinity
   double inf = numeric_limits<double>::infinity();
   int inf_int = numeric_limits<int>::infinity();

   //Initialize the routes
   vector<vector<vector<double>>> f(limit, vector<vector<double>>(len_values, vector<double> (len_N)));
   vector<vector<vector<double>>> phi(limit, vector<vector<double>>(len_values, vector<double> (len_N)));
   vector<vector<vector<int>>> p(limit, vector<vector<int>>(len_values, vector<int> (len_N)));
   vector<vector<vector<vector<int>>>> q_route(limit, vector<vector<vector<int>>>(len_values, vector<vector<int>> (len_N, vector<int>(0))));
   vector<vector<vector<vector<int>>>> q_route_2(limit, vector<vector<vector<int>>>(len_values, vector<vector<int>> (len_N, vector<int>(0))));




}

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
);
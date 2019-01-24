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
#include <list>
#include "lower_bounds.h"
#include "baldacci.h"

#include<cmath>
#include "prettyprint.hpp"
#include <bitset>



vector<DualSolution> construct_lower_bound_wrapper(
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
   vector<int> H,
   vector<int> capacities,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> geo_distance,
   vector<int> n_trucks,
   vector<vector<double>> penalties,
   double penalty_factor
);
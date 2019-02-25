
#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <random>
#include <functional>
#include<cmath>
#include "lower_bounds.h"
#include "baldacci.h"
#include "lower_bounds.cpp"
#include "VRPClass.cpp"

#include "baldacci.cpp"
template class std::set<tuple<int,int>>;
# include <chrono>
using  ms = chrono::milliseconds;
using get_time = chrono::steady_clock ;
#include "prettyprint.hpp"

#include "nlohmann/json.hpp"
#include <fstream>

// for convenience
using json = nlohmann::json;


int main(){

   VRP vrp = VRP("/Users/sergiocamelo/Dropbox/Sergio-Joann/Results/Jan242019-nopenalty/instances/daily/daily_cluster_780_day_12.json");





   int iterations_grad_m1 = 200;
   int iterations_grad_m2 = 100;
   int iterations_m2 = 3;
   double z_ub = 500000;
   int Delta = 2;
   int Delta_zero = Delta;
   int Delta_final = Delta;
   double gamma = 20000;
   double gamma_zero = - pow(10,-14)*z_ub;
   double gamma_final = 20000;
   double epsilon = 0.1;




   vector<DualSolution> lb = construct_lower_bound(iterations_grad_m1,
      iterations_grad_m2,
      iterations_m2,
      z_ub,
      Delta,
      Delta_zero,
      Delta_final,
      gamma,
      gamma_zero,
      gamma_final,
      epsilon,
      vrp
   );




   return 1;
}

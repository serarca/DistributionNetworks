
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

int main(int argc, char** argv){

   cout<<"Running instance: "<<argv[2]<<endl;

   string results_folder = argv[1];
   string instance_name = argv[2];
   string instances_folder = "instances/";

   VRP vrp = VRP_from_filename(results_folder+instances_folder+instance_name+".json", instance_name);
   vrp.folder = results_folder;
   cout<<"No. Farmers "<<vrp.len_N()<<endl;
   cout<<"No. Trucks "<<vrp.len_H()<<endl;



   int iterations_grad_m1 = 100;
   int iterations_grad_m2 = 150;
   int iterations_m2 = 3;
   double z_ub = 2000000;
   int Delta = 1000;
   int Delta_zero = Delta;
   int Delta_final = Delta;
   double gamma = 20000;
   double gamma_zero = - pow(10,-14)*z_ub;
   double gamma_final = 20000;
   double epsilon = 0.1;
   int limit = 15;


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
      vrp,
      limit
   );


   return 1;
}

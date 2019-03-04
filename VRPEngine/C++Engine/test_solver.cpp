
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
#include "utils/save_dual_solution.cpp"


// for convenience
using json = nlohmann::json;


int main(){

   string results = "/Users/sergiocamelo/Dropbox/Sergio-Joann/Results/Jan242019-nopenalty/";
   string folder = "instances/";
   string filename = "daily/daily_cluster_780_day_6";
   //filename = "daily/daily_cluster_780_day_12";
   VRP vrp = VRP_from_filename(results+folder+filename+".json", filename);
   vrp.folder = results;
   cout<<"No. Farmers "<<vrp.len_N()<<endl;


   int iterations_grad_m1 = 50;
   int iterations_grad_m2 = 100;
   int iterations_m2 = 3;
   double z_ub = 1000000;
   int Delta = 1000;
   int Delta_zero = Delta;
   int Delta_final = Delta;
   double gamma = 20000;
   double gamma_zero = - pow(10,-14)*z_ub;
   double gamma_final = 20000;
   double epsilon = 0.1;

   // //string iter = "1";
   // DualSolution partial_sol = read_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+".json");
   // // DualSolution partial_sol = read_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+"_partial.json");
   // partial_sol.initialize_routes(vrp.len_H());
   // DualSolution lb1 = optimize_lower_bound_M2(iterations_grad_m2, z_ub, Delta, Delta_zero, gamma, gamma_zero, epsilon, vrp, partial_sol);
   // save_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+"_lb1"+".json", lb1, vrp);
   //
   // DualSolution lb2 = optimize_lower_bound_M2(iterations_grad_m2, z_ub, Delta, Delta_zero, gamma, gamma_zero, epsilon, vrp, lb1);
   // save_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+"_lb2"+".json", lb2, vrp);



   // generateTruckMinRoutes(
   //    z_ub,
   //    16000,
   //    gamma_zero,
   //    80,
   //    partial_sol,
   //    vrp,
   //    - pow(10,-13) * z_ub,
   //    true
   // );

   // double truck_zero_gamma_guarantee = 0;
   // int Delta_zero_current = 1000;
   // bool termination = false;
   //
   // int h = 97;
   //list<SimpleRoute> routes = GENROUTE(z_ub, Delta_zero_current, gamma_zero, h, termination, truck_zero_gamma_guarantee, vrp, partial_sol);


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


   save_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+".json", lb[iterations_m2], vrp);


   // cout<<vrp.geo_distance[80][62]<<endl;
   // cout<<vrp.geo_distance[27][80]<<endl;
   // cout<<vrp.geo_distance[80][27]<<endl;
   //
   // cout<<vrp.geo_distance[62][80]<<endl;


   return 1;
}

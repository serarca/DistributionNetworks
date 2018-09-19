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
#include "baldacci.cpp"
#include "reduced_routes.cpp"

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

   std::ifstream i("/Users/sergiocamelo/Dropbox/Sergio-Joann/StandardizedData/instances/daily/daily_cluster_108_day_9.json");

   json j;
   i >> j;

   vector<int> N = j["N"];
   vector<int> H = j["H"];
   vector<int> quantities = j["quantities"];
   vector<int> capacities = j["capacities"];
   vector<int> n_trucks = j["n_trucks"];

   vector<vector<double>> geo_distance = j["distances"];
   cout<<geo_distance<<endl;

   double z_ub = 500000;


   int H_len = H.size();
   int N_len = N.size();

   cout<<"There is "<<H_len<<" trucks"<<endl;
   cout<<"There is "<<N_len<<" farmers"<<endl;
   cout<<"Capacities:"<<capacities<<endl;
   cout<<"Quantities:"<<quantities<<endl;




   vector<double> mu(H_len,0);
   vector<double> u(N_len,0);




   int iterations_grad_m1 = 200;
   int iterations_grad_m2 = 100;
   int iterations_m2 = 3;
   int Delta = 6000;
   int Delta_zero = Delta;
   int Delta_final = Delta;
   double gamma = 20000;
   double gamma_zero = - pow(10,-14)*z_ub;
   double gamma_final = 20000;
   double epsilon = 0.1;


   vector<list<SimpleRoute>> Routes = get_reduced_routes(
      z_ub,
      Delta,
      gamma,
      H,
      capacities,
      N,
      quantities,
      geo_distance,
      mu,
      u
   );


   print_sroute((Routes[0]).front());


   return 1;
}

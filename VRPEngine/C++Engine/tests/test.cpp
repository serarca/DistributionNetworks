
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

   int H_len = 1;
   int N_len = 10;
   vector<int> capacities(H_len, 10);
   vector<int> quantities(N_len, 1);
   vector<int> n_trucks(H_len, 1);


   vector<int> N(N_len,0);
   for (int i = 0; i<N_len; i++){
      N[i] = i;
   }
   vector<int> H(H_len,0);
   for (int i = 0; i<H_len; i++){
      H[i] = N_len+i;
   }
   vector<double> mu(H_len,0);
   vector<double> lamb(N_len,0);

   vector<double> x(N_len+H_len);
   vector<double> y(N_len+H_len);


   // First create an instance of an engine.
   random_device rnd_device;
   // Specify the engine and distribution.
   //mt19937 mersenne_engine(rnd_device());
   string str = "a";
   std::seed_seq seed1(str.begin(),str.end());
   mt19937 mersenne_engine(seed1);
   //mersenne_engine.seed(std::random_device()());
   ///mersenne_engine.seed(0);
   uniform_real_distribution<double> dist(0, 1);

   auto gen = std::bind(dist, mersenne_engine);
   generate(begin(x), end(x), gen);

   //mersenne_engine.seed(std::random_device()());
   str = "b";
   std::seed_seq seed2(str.begin(),str.end());
   mt19937 mersenne_engine_2(seed2);
   auto gen2 = std::bind(dist, mersenne_engine_2);
   generate(begin(y), end(y), gen2);

   // mu.push_back(0);
   // capacities.push_back(10);
   // n_trucks.push_back(1);
   // H.push_back(N_len+H_len);
   // x.push_back(x[N_len+H_len-1]);
   // y.push_back(x[N_len+H_len-1]);
   // H_len+=1;



   vector<vector<double>> geo_distance(N_len+H_len, vector<double>(N_len+H_len));
   for (int i = 0; i < N_len+H_len; i++){
      for (int j = 0; j < N_len+H_len; j++){
         geo_distance[i][j] = sqrt(pow(x[i]-x[j],2)+pow(y[i]-y[j],2));
      }
   }

   //LowerBound lb = lower_bound_(H, capacities, N, quantities, distance_dict, mu, lamb);
   //cout<<lb.z_lb<<endl;

   /*
   int Delta = 700;
   double gamma = 12;
   int h = H[0];
   int capacity = capacities[0];
   string direction = "left";
   //vector<list<Path>> p = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, direction);
   auto start = get_time::now();
   vector<list<SimpleRoute>> r = GENROUTE(Delta, gamma, h, capacity, N, quantities, distance_dict);
   auto end = get_time::now();
   auto diff = end - start;
   cout<<"Elapsed time is :  "<< chrono::duration_cast<ms>(diff).count()<<" ms "<<endl;
   print_SRoutes(r);
   */
   /*
   // Debugging first lower bounds
   double z_ub = 20;
   int iterations = 100;
   double epsilon = 0.1;
   DualSolution sol = lower_bound_optimizer_M1(iterations, z_ub, epsilon, H, capacities, N, quantities,geo_distance);
   cout<<"Bound Christofides:"<<sol.z_lb<<endl;
   sol.routes.clear();


   // Debugging the lower_bound_M2
   int Delta = 1000;
   int Delta_zero = 1000;
   double gamma = 5;
   double gamma_zero = 5;
   int len_N = N.size();
   int len_H = H.size();
   int sub_iterations = 100;
   int Delta_final = 1000;
   double gamma_final = 5;

   int iter_m2 = 5;
   for (int iter_2 = 0; iter_2<iter_m2; iter_2++){
      sol = optimize_lower_bound_M2(sub_iterations, z_ub, Delta, Delta_zero, gamma, gamma_zero, epsilon, H, capacities, N, quantities, geo_distance, sol.v, sol.lamb, sol.u, sol.routes);
      cout<<"Bound M2:"<<sol.z_lb<<endl;
   }

   // Define the vector of routes for each truck
   vector<list<SimpleRoute>> FinalRoutes(len_H);

   // Calculate reduced costs
   vector<vector<double>> distance_dict = reduced_cost_matrix(geo_distance, sol.u, sol.v);

   // We start by generating routes for all of the trucks
   for (auto h:H){
      FinalRoutes[h - len_N] = GENROUTE(Delta_final, gamma_final, h, capacities[h - len_N], N, quantities, distance_dict, geo_distance);
   }

   // Now we add routes that come from a previous iteration of the algorithm
   if (sol.routes.size()>0){
      for (int i = 0; i< len_H; i++){
         for (auto old_route:(sol.routes)[i]){
            bool add = true;
            for (auto new_route:FinalRoutes[i]){
               if (new_route.path == old_route.path){
                  add = false;
                  break;
               }
            }
            if (add){
               FinalRoutes[i].push_back(old_route);
            }
         }
      }
   }
   cout<<"Routes per truck"<<endl;
   for (int i = 0; i< len_H; i++){
      cout<<i+len_H<<":"<<FinalRoutes[i].size()<<endl;
   }
   */
   int iterations_grad_m1 = 200;
   int iterations_grad_m2 = 100;
   int iterations_m2 = 1;
   double z_ub = 20;
   int Delta = 3000;
   int Delta_zero = 1000;
   int Delta_final = 3000;
   double gamma = 1;
   double gamma_zero = - pow(10,-14);
   double gamma_final = 1;
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
      H,
      capacities,
      N,
      quantities,
      geo_distance,
      n_trucks
   );
   cout<<lb[0].routes.size()<<endl;


   return 1;
}

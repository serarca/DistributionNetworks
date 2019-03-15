
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
#include "q_paths.cpp"


// for convenience
using json = nlohmann::json;


int main(){

   string results = "/Users/sergiocamelo/Dropbox/Sergio-Joann/Results/Jan242019-nopenalty/";
   string folder = "instances/";
   string filename = "spatial/spatial_day_11";
   //filename = "daily/daily_cluster_780_day_12";
   cout<<results+folder+filename+".json"<<endl;
   VRP vrp = VRP_from_filename(results+folder+filename+".json", filename);
   vrp.folder = results;
   cout<<"No. Farmers "<<vrp.len_N()<<endl;

   cout<<vrp.folder+"dual_solutions/"+vrp.name+".json"<<endl;
   DualSolution partial_sol = read_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+"_iter_0.json");
   partial_sol.calc_reduced_distances(vrp.geo_distance);

   auto distances = partial_sol.reduced_distances;

   string direction = "left";

   int limit = 15;
   int h = 70;
   QPaths qpaths = construct_q_paths_(
      vrp,
      h,
      distances,
      direction
   );
   //cout<<qpaths.q_route[80][10]<<endl;

   QPathsLim qpaths_lim = construct_q_paths_lim_(
      vrp,
      h,
      distances,
      direction,
      limit
   );

   // Check entries
   // cout<<qpaths_lim.f[0][80][10]<<endl;
   // cout<<qpaths_lim.q_route[0][80][10]<<endl;
   // cout<<qpaths_lim.f[1][80][10]<<endl;
   // cout<<qpaths_lim.q_route[1][80][10]<<endl;
   // cout<<qpaths_lim.f[2][80][10]<<endl;
   // cout<<qpaths_lim.q_route[2][80][10]<<endl;
   // cout<<qpaths_lim.f[3][80][10]<<endl;
   // cout<<qpaths_lim.q_route[3][80][10]<<endl;
   // cout<<qpaths_lim.f[4][80][10]<<endl;
   // cout<<qpaths_lim.q_route[4][80][10]<<endl;
   // cout<<qpaths_lim.f[5][80][10]<<endl;
   // cout<<qpaths_lim.q_route[5][80][10]<<endl;
   // cout<<qpaths_lim.f[6][80][10]<<endl;
   // cout<<qpaths_lim.q_route[6][80][10]<<endl;


   double inf = numeric_limits<double>::infinity();
   int inf_int = numeric_limits<int>::max();
   // Verify the length of the paths
   bool pass = true;
   for (int lim = 0; lim < limit; lim++){
      for (int l = 0; l < vrp.possible_v[h-vrp.len_N()].values.size(); l++){
         for (int end = 0; end < vrp.len_N(); end++){
            // Verify the length of the paths
            if (qpaths_lim.q_route[lim][l][end].size() == 0 && qpaths_lim.f[lim][l][end] != inf){
               cout<<"Error1"<<endl;
            }
            if (qpaths_lim.q_route[lim][l][end].size() != 0 && qpaths_lim.f[lim][l][end] == inf){
               cout<<"Error2"<<endl;
            }
            if (qpaths_lim.q_route_2[lim][l][end].size() == 0 && qpaths_lim.phi[lim][l][end] != inf){
               cout<<"Error3"<<endl;
            }
            if (qpaths_lim.q_route_2[lim][l][end].size() != 0 && qpaths_lim.phi[lim][l][end] == inf){
               cout<<"Error4"<<endl;
            }
            if (qpaths_lim.q_route_2[lim][l][end].size() != 0 && qpaths_lim.q_route_2[lim][l][end].size() != lim+2){
               cout<<"Error5"<<endl;
            }
            if (qpaths_lim.q_route[lim][l][end].size() != 0 && qpaths_lim.q_route[lim][l][end].size() != lim+2){
               cout<<"Error6"<<endl;
            }
            // Correct length and no repetition of vertices
            double length = 0.0;
            int load = 0;
            if (qpaths_lim.q_route[lim][l][end].size()!=0){
               for (int i = 0; i<qpaths_lim.q_route[lim][l][end].size()-1; i++){
                  if (direction == "left"){
                     length += distances[qpaths_lim.q_route[lim][l][end][i]][qpaths_lim.q_route[lim][l][end][i+1]];
                  }else{
                     length += distances[qpaths_lim.q_route[lim][l][end][i+1]][qpaths_lim.q_route[lim][l][end][i]];
                  }                  if (i<qpaths_lim.q_route[lim][l][end].size()-2){
                     if (qpaths_lim.q_route[lim][l][end][i]==qpaths_lim.q_route[lim][l][end][i+2] ||
                     qpaths_lim.q_route[lim][l][end][i]==qpaths_lim.q_route[lim][l][end][i+1]){
                        cout<<"Error11"<<endl;
                     }
                  }
               }
               if (abs(length-qpaths_lim.f[lim][l][end])> 1.0e-10){
                  cout<<"Error10"<<endl;
               }
               for (int i = 0; i<qpaths_lim.q_route[lim][l][end].size(); i++){
                  if (qpaths_lim.q_route[lim][l][end][i]<vrp.len_N()){
                     load += vrp.quantities[qpaths_lim.q_route[lim][l][end][i]];
                  }
               }
               if (load!=vrp.possible_v[h-vrp.len_N()].values[l]){
                  cout<<"Error_Quant"<<endl;
               }
            }
            // Correct length and no repetition of vertices
            load = 0;
            length = 0.0;
            if (qpaths_lim.q_route_2[lim][l][end].size()!=0){
               for (int i = 0; i<qpaths_lim.q_route_2[lim][l][end].size()-1; i++){
                  if (direction == "left"){
                     length += distances[qpaths_lim.q_route_2[lim][l][end][i]][qpaths_lim.q_route_2[lim][l][end][i+1]];
                  }else{
                     length += distances[qpaths_lim.q_route_2[lim][l][end][i+1]][qpaths_lim.q_route_2[lim][l][end][i]];
                  }
                  if (i<qpaths_lim.q_route_2[lim][l][end].size()-2){
                     if (qpaths_lim.q_route_2[lim][l][end][i]==qpaths_lim.q_route_2[lim][l][end][i+2] ||
                     qpaths_lim.q_route_2[lim][l][end][i]==qpaths_lim.q_route_2[lim][l][end][i+1]){
                        cout<<"Error11"<<endl;
                     }
                  }
               }
               if (abs(length-qpaths_lim.phi[lim][l][end])> 1.0e-10){
                  cout<<"Error10"<<endl;
               }
               for (int i = 0; i<qpaths_lim.q_route[lim][l][end].size(); i++){
                  if (qpaths_lim.q_route[lim][l][end][i]<vrp.len_N()){
                     load += vrp.quantities[qpaths_lim.q_route[lim][l][end][i]];
                  }
               }
               if (load!=vrp.possible_v[h-vrp.len_N()].values[l]){
                  cout<<"Error_Quant"<<endl;
                  cout<<load<<endl;
                  cout<<vrp.possible_v[h-vrp.len_N()].values[l]<<endl;
               }
            }

         }
      }
   }

   //q_route agrees with q_route_lim
   for (int l = 0; l < vrp.possible_v[h-vrp.len_N()].values.size(); l++){
      for (int end = 0; end < vrp.len_N(); end++){
         int lim = qpaths.q_route[l][end].size();
         if (lim > 0 && lim-2 < limit){
            // cout<<endl;
            // cout<<"lim:"<<lim<<endl;
            // cout<<"qpaths_lim.f(len):"<<qpaths_lim.f.size()<<endl;
            // cout<<"l:"<<l<<endl;
            // cout<<"qpaths(l):"<<qpaths.f.size()<<endl;
            // cout<<"qpaths_lim.f[lim-2](l):"<<qpaths_lim.f[lim-2].size()<<endl;
            // cout<<"end:"<<end<<endl;
            // cout<<"qpaths.f[l](end):"<<qpaths.f[l].size()<<endl;
            // cout<<"qpaths_lim.f[lim-2][l](end):"<<qpaths_lim.f[lim-2][l].size()<<endl;
            if (qpaths.f[l][end] != qpaths_lim.f[lim-2][l][end]){
               cout<<"Error7"<<endl;
            }
            if (qpaths.q_route[l][end] != qpaths_lim.q_route[lim-2][l][end]){
               cout<<"Error8"<<endl;
            }
         }
      }
   }

   for (int l = 0; l < vrp.possible_v[h-vrp.len_N()].values.size(); l++){
      for (int end = 0; end < vrp.len_N(); end++){
         // Verify the length of the paths
         if (qpaths.q_route[l][end].size() == 0 && qpaths.f[l][end] != inf){
            cout<<"Error1"<<endl;
         }
         if (qpaths.q_route[l][end].size() != 0 && qpaths.f[l][end] == inf){
            cout<<"Error2"<<endl;
         }
         if (qpaths.q_route_2[l][end].size() == 0 && qpaths.phi[l][end] != inf){
            cout<<"Error3"<<endl;
         }
         if (qpaths.q_route_2[l][end].size() != 0 && qpaths.phi[l][end] == inf){
            cout<<"Error4"<<endl;
         }

         // Correct length and no repetition of vertices
         double length = 0.0;
         if (qpaths.q_route[l][end].size()!=0){
            for (int i = 0; i<qpaths.q_route[l][end].size()-1; i++){
               if (direction == "left"){
                  length += distances[qpaths.q_route[l][end][i]][qpaths.q_route[l][end][i+1]];
               }else{
                  length += distances[qpaths.q_route[l][end][i+1]][qpaths.q_route[l][end][i]];
               }                  if (i<qpaths.q_route[l][end].size()-2){
                  if (qpaths.q_route[l][end][i]==qpaths.q_route[l][end][i+2] ||
                  qpaths.q_route[l][end][i]==qpaths.q_route[l][end][i+1]){
                     cout<<"Error11"<<endl;
                  }
               }
            }
            if (abs(length-qpaths.f[l][end])> 1.0e-10){
               cout<<"Error10"<<endl;
            }
         }
         // Correct length and no repetition of vertices

         length = 0.0;
         if (qpaths.q_route_2[l][end].size()!=0){
            for (int i = 0; i<qpaths.q_route_2[l][end].size()-1; i++){
               if (direction == "left"){
                  length += distances[qpaths.q_route_2[l][end][i]][qpaths.q_route_2[l][end][i+1]];
               }else{
                  length += distances[qpaths.q_route_2[l][end][i+1]][qpaths.q_route_2[l][end][i]];
               }
               if (i<qpaths.q_route_2[l][end].size()-2){
                  if (qpaths.q_route_2[l][end][i]==qpaths.q_route_2[l][end][i+2] ||
                  qpaths.q_route_2[l][end][i]==qpaths.q_route_2[l][end][i+1]){
                     cout<<"Error11"<<endl;
                  }
               }
            }
            if (abs(length-qpaths.phi[l][end])> 1.0e-10){
               cout<<"Error10"<<endl;
            }
         }

      }
   }

   //Now we calculate the q_routes
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();
   QRoutes qroutes = construct_q_routes_(vrp, h, distances);
   end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end - start;

   std::cout << elapsed_seconds.count() << "s\n";
   start = std::chrono::system_clock::now();
   QRoutes qroutes_lim = construct_q_routes_lim_(vrp, h, distances, limit);
   end = std::chrono::system_clock::now();
   elapsed_seconds = end - start;
   std::cout << elapsed_seconds.count() << "s\n";

   //Debug stuff

   cout<<qroutes.psi_route[80][10]<<endl;
   cout<<qroutes_lim.psi_route[80][10]<<endl;
   cout<<qroutes.psi_route[80][12]<<endl;
   cout<<qroutes_lim.psi_route[80][12]<<endl;

   //Verify that everything is fine
   for (int l = 0; l < vrp.possible_v[h-vrp.len_N()].values.size(); l++){
      for (int end = 0; end < vrp.len_N(); end++){
         // Lim routes have the required length
         auto q_route = qroutes.psi_route[l][end];
         auto q_route_lim = qroutes_lim.psi_route[l][end];
         if (q_route_lim.size()>0){
            if (q_route_lim.size()>limit+2){
               cout<<"Error:length"<<endl;
            }
         }
         if (qroutes.psi[l][end]>qroutes_lim.psi[l][end]){
            cout<<"Error:value"<<endl;
            cout<<(qroutes.psi[l][end]-qroutes_lim.psi[l][end])<<endl;
         }
         if (qroutes.psi[l][end]!=qroutes_lim.psi[l][end]){
            cout<<endl;
            cout<<qroutes.psi_route[l][end]<<endl;
            cout<<qroutes_lim.psi_route[l][end]<<endl;
            cout<<(qroutes.psi[l][end]-qroutes_lim.psi[l][end])<<endl;
         }
      }
   }

   // Test the lower bounds
   double epsilon = 0.1;
   int iterations = 10;
   double z_ub = 1000000;

   DualSolution sol1 = lower_bound_optimizer_M1(
      iterations,
      z_ub,
      epsilon,
      vrp,
      inf_int
   );

   DualSolution sol2 = lower_bound_optimizer_M1(
      iterations,
      z_ub,
      epsilon,
      vrp,
      15
   );




   return 1;
}

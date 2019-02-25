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


#include "nlohmann/json.hpp"
#include <fstream>
// for convenience
using json = nlohmann::json;


VRP read_VRP(std::string filename, double penalty_factor = 0.0){

   std::ifstream i("/Users/sergiocamelo/Dropbox/Sergio-Joann/Results/2018-10-20_17:58/instances/daily/daily_cluster_591_day_11.json");

   json j;
   i >> j;

   vector<int> N = j["N"];
   vector<int> H = j["H"];
   vector<int> quantities = j["quantities"];
   vector<int> capacities = j["capacities"];
   vector<int> n_trucks = j["n_trucks"];

   vector<vector<double>> geo_distance = j["distances"];

   // Extract penalties
   auto it_penalties = j.find("penalties");
   vector<vector<double>> penalties;
   if (it_penalties != j.end()) {
      vector<vector<double>> penalties_extracted = *it_penalties;
      penalties = penalties_extracted;
   }

   int H_len = H.size();
   int N_len = N.size();

   cout<<"There is "<<H_len<<" trucks"<<endl;
   cout<<"There is "<<N_len<<" farmers"<<endl;
   cout<<"Capacities:"<<capacities<<endl;
   cout<<"Quantities:"<<quantities<<endl;

   VRP vrp = VRP();
   if (penalties.size() > 0) {
      vrp = VRP(H, capacities, N, quantities, geo_distance, n_trucks, penalties);
   } else {
      vrp = VRP(H, capacities, N, quantities, geo_distance, n_trucks);
   }
   vrp.penalty_factor = penalty_factor;


   return vrp;
}

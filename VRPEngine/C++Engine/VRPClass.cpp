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


#include<cmath>
#include "prettyprint.hpp"
#include <bitset>


#include "nlohmann/json.hpp"
#include <fstream>
// for convenience
using json = nlohmann::json;

using namespace std;


class VRP {

public:

   vector<int> H;

   vector<int> capacities;
   vector<int> N;
   vector<int> quantities;
   vector<vector<double>> geo_distance;
   vector<int> n_trucks;
   vector<vector<double>> penalties;
   bool penalized;
   double penalty_factor;
   string name;
   string folder;

   VRP(){}


   VRP(
      vector<int> H_,
      vector<int> capacities_,
      vector<int> N_,
      vector<int> quantities_,
      vector<vector<double>> geo_distance_,
      vector<int> n_trucks_,
      vector<vector<double>> penalties_,
      double penalty_factor_ = 0.0
   ){
      H = H_;
      capacities = capacities_;
      N = N_;
      quantities = quantities_;
      geo_distance = geo_distance_;
      n_trucks = n_trucks_;
      penalties = penalties_;
      penalized = true;
      penalty_factor = penalty_factor_;
   }

   VRP(
      vector<int> H_,
      vector<int> capacities_,
      vector<int> N_,
      vector<int> quantities_,
      vector<vector<double>> geo_distance_,
      vector<int> n_trucks_
   ){
      H = H_;
      capacities = capacities_;
      N = N_;
      quantities = quantities_;
      geo_distance = geo_distance_;
      n_trucks = n_trucks_;
      penalized = false;
      penalty_factor = 0.0;
   }

   int len_H(){
      return H.size();
   }
   int len_N(){
      return N.size();
   }

   vector<int> V(){
      vector<int> V;
      V.reserve( N.size() + H.size() ); // preallocate memory
      V.insert( V.end(), N.begin(), N.end() );
      V.insert( V.end(), H.begin(), H.end() );
      return V;
   }

};

VRP VRP_from_filename(string filename, string name_="", double penalty_factor_ = 0.0){
   std::ifstream i(filename);
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

   VRP vrp;

   if (penalties.size() > 0) {
      vrp = VRP(H,capacities,N,quantities,geo_distance,n_trucks, penalties, penalty_factor_);
   } else {
      vrp = VRP(H,capacities,N,quantities,geo_distance,n_trucks);
   }
   vrp.name = name_;

   return vrp;



}

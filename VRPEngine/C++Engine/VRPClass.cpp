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

// A struct with posible values and their inverse maps
struct PossibleValues {
   vector<int> values;
   map<int,int> values_pos;
};

PossibleValues possible_values(vector<int>& quantities, int truck_capacity){
   vector<int> values(truck_capacity);

   for (int i=0; i< truck_capacity; i++){
      values[i] = i + 1;
   }
   map<int, int> values_pos;
   for (int i = 0; i < (int) values.size(); i++){
      values_pos[values[i]] = i;
   }
   PossibleValues possible;
   possible.values = values;
   possible.values_pos = values_pos;
   return possible;
}

vector<vector<double>> penalize_distance(
   vector<vector<double>> distance_dict,
   vector<double> penalties,
   int len_N,
   int len_H,
   double penalty_factor
){
   for (int n = 0; n < len_N; n++){
      for (int v = 0; v < len_N + len_H; v++){
         if (v != n){
            distance_dict[n][v] += penalties[n]/2.0 * penalty_factor;
            distance_dict[v][n] += penalties[n]/2.0 * penalty_factor;
         }
      }
   }
   return distance_dict;
}



class VRP {

public:

   vector<int> H;

   vector<int> capacities;
   vector<int> N;
   vector<int> quantities;
   vector<vector<double>> geo_distance;
   vector<vector<vector<double>>> truck_distances;
   vector<int> n_trucks;
   vector<vector<double>> penalties;
   bool penalized;
   double penalty_factor;
   string name;
   string folder;
   vector<PossibleValues> possible_v;
   string mapping;

   VRP(){}


   VRP(
      vector<int> H_,
      vector<int> capacities_,
      vector<int> N_,
      vector<int> quantities_,
      vector<vector<vector<double>>> truck_distances_,
      vector<int> n_trucks_
   ){
      // Set up basic quantities
      H = H_;
      capacities = capacities_;
      N = N_;
      quantities = quantities_;
      truck_distances = truck_distances_;
      n_trucks = n_trucks_;

      // Initialize possible values
      initialize_possible_values();

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

   void initialize_possible_values(){
      possible_v = vector<PossibleValues>(len_H());
      for(int i = 0; i < len_H(); i++){
         possible_v[i] = possible_values(quantities, capacities[i]);
      }
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
   string mapping = j["mapping"].dump();

   vector<vector<double>> geo_distance = j["distances"];

   // Extract penalties
   auto it_penalties = j.find("penalties");
   vector<vector<double>> penalties;
   if (it_penalties != j.end()) {
      vector<vector<double>> penalties_extracted = *it_penalties;
      penalties = penalties_extracted;
   }

   VRP vrp;

   // Construct truck_distances
   auto it_truck_distances = j.find("truck_distances");
   vector<vector<vector<double>>> truck_distances;
   if (it_truck_distances != j.end()) {
      vector<vector<vector<double>>> distances_extracted = *it_truck_distances;
      truck_distances = distances_extracted;
   } else {
      //Construct truck distances using geo_distances
      for(int i=0; i<H.size(); i++){
         truck_distances.push_back(geo_distance);
      }
   }

   vrp = VRP(H,capacities,N,quantities,truck_distances,n_trucks);

   if (penalties.size() > 0 && penalty_factor_!=0.0) {
      vrp.penalized = true;
      for(int i=0; i<vrp.len_H(); i++){
         vrp.truck_distances[i] = penalize_distance(vrp.truck_distances[i], penalties[i], N.size(), H.size(), penalty_factor_);
      }
      vrp.penalties = penalties;
      vrp.penalty_factor = penalty_factor_;
   } else {
      vrp.penalized = false;
   }

   // Add other quantities
   vrp.geo_distance = geo_distance;
   vrp.name = name_;
   vrp.mapping = mapping;

   return vrp;



}

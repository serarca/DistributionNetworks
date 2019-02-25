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


   VRP(string filename, double penalty_factor_ = 0.0, string name_=""){
      name = name_;

      std::ifstream i(filename);

      json j;
      i >> j;

      N = j["N"].get<vector<int>>();
      H = j["H"].get<vector<int>>();
      quantities = j["quantities"].get<vector<int>>();
      capacities = j["capacities"].get<vector<int>>();
      n_trucks = j["n_trucks"].get<vector<int>>();

      geo_distance = j["distances"].get<vector<vector<double>>>();

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


      if (penalties.size() > 0) {
         penalized = true;
      } else {
         penalized = false;
      }
      penalty_factor = penalty_factor_;

      cout<< H;

      cout<< capacities<<endl;
      cout<< N<<endl;
      cout<< quantities<<endl;
      //cout<< geo_distance<<endl<<endl;
      cout<< n_trucks<<endl;
      cout<< penalties<<endl;
      cout<< penalized<<endl;
      cout<< penalty_factor<<endl;
      cout<< name<<endl;


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

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
#include "VRPClass.cpp"
#include "lower_bounds.h"


// for convenience
using json = nlohmann::json;

using namespace std;

class PrimalSolution{

public:
   VRP vrp;
   vector<SimpleRoute> routes;
   double z_lb;
   int load;

   PrimalSolution(VRP& vrp_){
      vrp = vrp_;
      vector<SimpleRoute> routes_;
      for (int i = 0; i < vrp.len_H(); i++){
         SimpleRoute route;
         routes_.push_back(route);
      }
      routes = routes_;
   }

   void update_routes(){
      for (int i = 0; i < vrp.len_H(); i++){
         SimpleRoute& route = routes[i];
         route.geo_cost = 0;
         route.load = 0;
         // create path
         vector<int> v_path{ std::begin(route.path), std::end(route.path) };
         if (v_path.size()>0){
            for (int j=0; j<v_path.size()-1;j++){
               route.geo_cost += vrp.geo_distance[v_path[j]][v_path[j+1]];
            }
            for (int j=0; j<v_path.size();j++){
               if (v_path[j] < vrp.len_N()){
                  route.load += vrp.quantities[v_path[j]];
               }
            }
         }
         route.cost = route.geo_cost;
      }
   }

   bool verify_solution(){
      z_lb = 0;
      load = 0;
      vector<int> served_farmers(vrp.len_N(),0);
      for (int i = 0; i < vrp.len_H(); i++){
         SimpleRoute& route = routes[i];
         z_lb += route.cost;
         load += route.load;
         if (vrp.capacities[i]<route.load){
            cout<<"Violating loads"<<endl;
            return false;
         }

         vector<int> v_path{ std::begin(route.path), std::end(route.path) };
         if (v_path.size()>0){
            for (int j=0; j<v_path.size();j++){
               if (v_path[j] < vrp.len_N()){
                  served_farmers[v_path[j]] += 1;
               }
            }
         }
      }
      for (int i = 0; i<served_farmers.size(); i++){
         if (served_farmers[i]!= 1){
            cout<<"Violating farmers"<<endl;
            return false;
         }
      }
      cout<<"Total Cost: "<<z_lb<<endl;
      cout<<"Total Load: "<<load<<endl;
      return true;
   }

   void save_solution(string filename){
      json j;

      update_routes();
      if (!verify_solution()){
         j["status"] = "verified";
      } else {
         j["z_lb"] = z_lb;

         for (int i = 0; i < vrp.len_H(); i++){
            SimpleRoute& route = routes[i];
            vector<int> v_path{ std::begin(route.path), std::end(route.path) };
            j["routes"].push_back(v_path);
         }

      }

      std::ofstream o(filename);
      o << std::setw(4) << j << std::endl;



   }



};

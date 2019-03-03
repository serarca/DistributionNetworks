#pragma once


#include "../lower_bounds.h"

#include "../nlohmann/json.hpp"
#include <fstream>

// Saves dual solution to json file
void save_dual_solution(string file_name, DualSolution sol){
   json j;
   j["z_lb"] = sol.z_lb;
   j["lamb"] = sol.lamb;
   j["u"] =  sol.u;
   j["v"] = sol.v;

   cout<<file_name<<endl;

   std::ofstream file(file_name);
   file << j;
}

// Saves dual solution to json file
DualSolution read_dual_solution(string file_name){

   DualSolution sol;

   std::ifstream i(file_name);
   json j;
   i >> j;

   double z_lb = j["z_lb"];
   vector<double> lamb = j["lamb"];
   vector<double> u = j["u"];
   vector<double> v = j["v"];

   sol.z_lb = z_lb;
   sol.lamb = lamb;
   sol.u = u;
   sol.v = v;

   return sol;
}

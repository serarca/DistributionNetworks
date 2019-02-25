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

   std::ofstream file(file_name);
   file << j;
}

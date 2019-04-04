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
#include "utils/save_dual_solution.cpp"
#include "primal_solver.cpp"
#include "ParametersClass.cpp"



// for convenience
using json = nlohmann::json;

int main(int argc, char** argv){

   // Take time
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   cout<<"Running instance: "<<argv[2]<<endl;

   RouteParameters parameters;


   string results_folder = argv[1];
   string instance_name = argv[2];
   parameters.Delta = atoi(argv[3]);
   parameters.gamma = stod(argv[4]);
   parameters.z_ub = 1000000;
   parameters.route_limit = 15;

   cout<<"Running with Delta: "<<parameters.Delta<<endl;
   cout<<"Running with Gamma: "<<parameters.gamma<<endl;




   string instances_folder = "instances/";

   VRP vrp = VRP_from_filename(results_folder+instances_folder+instance_name+".json", instance_name);
   vrp.folder = results_folder;

   cout<<"No. Farmers "<<vrp.len_N()<<endl;
   cout<<"No. Trucks "<<vrp.len_H()<<endl;


   cout<< vrp.folder+"dual_solutions/"+vrp.name+"_iter_3.json" <<endl;
   string file_write = vrp.folder+"bash_c/spatial_primal/"+vrp.name+"_"+to_string(parameters.Delta);
   cout<<file_write<<endl;
   DualSolution sol = read_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+"_iter_3.json");
   PrimalSolution solution = primal_solution(vrp, sol, parameters);

   solution.save_solution(file_write);


   end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end - start;
   std::time_t start_time = std::chrono::system_clock::to_time_t(start);
   std::time_t end_time = std::chrono::system_clock::to_time_t(end);

   std::cout << "Started computation at " << std::ctime(&start_time)
   << "Finished computation at " << std::ctime(&end_time)
   << "Elapsed time: " << elapsed_seconds.count() << "s\n";


   return 1;
}

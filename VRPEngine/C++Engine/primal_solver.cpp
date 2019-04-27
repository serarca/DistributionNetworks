#pragma once

#include <iostream>
#include <fstream>
#include <gurobi_c++.h>
using namespace std;
#include "VRPClass.cpp"
#include "lower_bounds.cpp"
#include "baldacci.cpp"
#include "ParametersClass.cpp"
#include "PrimalClass.cpp"


PrimalSolution primal_solution(VRP& vrp, DualSolution& ds, RouteParameters& parameters){


   double inf = numeric_limits<double>::infinity();

   // Construct the routes from the dual solution
   cout<<"Updating Routes"<<endl;
   update_routes(
      ds,
      vrp,
      parameters.Delta,
      parameters.gamma,
      parameters.z_ub,
      inf,
      false,
      parameters.route_limit
   );
   cout<<"Finished updating Routes"<<endl;

   PrimalSolution solution(vrp);


   try {

   GRBEnv env = GRBEnv();
   GRBModel model = GRBModel(env);

   // Create variables
   vector<vector<GRBVar>> vars(vrp.len_H(), vector<GRBVar>(0));
   // Create farmer expressions
   vector<GRBLinExpr> farmer_expr(vrp.len_N(),GRBLinExpr(0.0));
   vector<GRBLinExpr> truck_expr(vrp.len_H(),GRBLinExpr(0.0));

   // Go through all routes
   for (int i = 0; i < ds.routes.size(); i++){
      int j = 0;
      for (auto route:ds.routes[i]){
         GRBVar var = model.addVar(0.0, 1.0, route.truck_cost, GRB_BINARY, to_string(i)+"_"+to_string(j));
         vars[i].push_back(var);
         for (auto farmer:route.path){
            if (farmer<vrp.len_N()){
               farmer_expr[farmer] += var;
            }
         }
         truck_expr[i] += var;
         j++;
      }
   }

   // Add Constraints
   for (int i = 0; i < vrp.len_N(); i++){
      model.addConstr(farmer_expr[i] == 1);
   }
   for (int i = 0; i < vrp.len_H(); i++){
      model.addConstr(truck_expr[i] <= 1);
   }


   // Set time limit
   model.getEnv().set(GRB_DoubleParam_TimeLimit, 7000);

   // Optimize model
   model.optimize();


   // Get status
   int optimstatus = model.get(GRB_IntAttr_Status);

   if (optimstatus == GRB_INFEASIBLE) {
      solution.feasible = false;
      cout << "Model is infeasible" << endl;
   }else{
      solution.feasible = true;
   }

   cout << "UB: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
   cout << "LB: " << ds.z_lb << endl;
   cout << "Gap: "<< (model.get(GRB_DoubleAttr_ObjVal) - ds.z_lb)/model.get(GRB_DoubleAttr_ObjVal)*100<<endl;


   // Extract solution
   for (int i = 0; i < ds.routes.size(); i++){
      int j = 0;
      for (auto route:ds.routes[i]){
         int value = round(vars[i][j].get(GRB_DoubleAttr_X));
         if (value == 1){
            solution.routes[i] = route;
         }
         j++;
      }
   }
   solution.z_lb = ds.z_lb;



   } catch(GRBException e) {
   cout << "Error code = " << e.getErrorCode() << endl;
   cout << e.getMessage() << endl;
   } catch(...) {
   cout << "Exception during optimization" << endl;
   }

   return solution;

}

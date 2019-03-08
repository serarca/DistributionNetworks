#pragma once


#include <gurobi_c++.h>
using namespace std;
#include "VRPClass.cpp"
#include "lower_bounds.cpp"
#include "baldacci.cpp"
#include "ParametersClass.cpp"

void primal_solution(VRP& vrp, DualSolution& ds, Parameters& parameters){


   double inf = numeric_limits<double>::infinity();

   // Construct the routes from the dual solution
   update_routes(
      ds,
      vrp,
      parameters.Delta_final,
      parameters.gamma_final,
      parameters.z_ub,
      inf,
      false
   );


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
         GRBVar var = model.addVar(0.0, 1.0, route.geo_cost, GRB_CONTINUOUS, to_string(i)+"_"+to_string(j));
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


   // Optimize model
   model.optimize();

   // cout << x.get(GRB_StringAttr_VarName) << " "
   // << x.get(GRB_DoubleAttr_X) << endl;
   // cout << y.get(GRB_StringAttr_VarName) << " "
   // << y.get(GRB_DoubleAttr_X) << endl;
   // cout << z.get(GRB_StringAttr_VarName) << " "
   // << z.get(GRB_DoubleAttr_X) << endl;

   cout << "UB: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
   cout << "LB: " << ds.z_lb << endl;
   cout << "Gap: "<< (model.get(GRB_DoubleAttr_ObjVal) - ds.z_lb)/model.get(GRB_DoubleAttr_ObjVal)*100<<endl;

   } catch(GRBException e) {
   cout << "Error code = " << e.getErrorCode() << endl;
   cout << e.getMessage() << endl;
   } catch(...) {
   cout << "Exception during optimization" << endl;
   }

}

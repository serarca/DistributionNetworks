#include <gurobi_c++.h>
using namespace std;
#include "VRPClass.cpp"


void create_cuts(VRP vrp){

   try {
   GRBEnv env = GRBEnv();

   GRBModel model = GRBModel(env);

   // Create variables
   vector<vector<GRBVar>> vars(vrp.len_N() + vrp.len_H(), vector<GRBVar>(vrp.len_N() + vrp.len_H()));

   for ( auto &i : vrp.V()){
      for ( auto &j : vrp.V()){
         // Check at least one node is a farmer
         if (i<vrp.len_N() || j<vrp.len_N()){
            double upper_bound;
            vars[i][j] = model.addVar(0.0, 1.0, vrp.geo_distance[i][j], GRB_CONTINUOUS, to_string(i)+"_"+to_string(j));
         }
      }
   }

   for ( auto &i : vrp.N){
      GRBLinExpr constr_ent = 0.0;
      GRBLinExpr constr_exi = 0.0;
      for ( auto &j : vrp.V()){
         constr_ent += vars[j][i];
         constr_exi += vars[i][j];
      }
      model.addConstr(constr_ent == 1.0);
      model.addConstr(constr_exi == 1.0);
   }

   GRBVar x = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x");
   GRBVar y = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y");
   GRBVar z = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "z");

   // Set objective: maximize x + y + 2 z

   model.setObjective(x + y + 2 * z, GRB_MAXIMIZE);

   // Add constraint: x + 2 y + 3 z <= 4

   model.addConstr(x + 2 * y + 3 * z <= 4, "c0");

   // Add constraint: x + y >= 1

   model.addConstr(x + y >= 1, "c1");

   // Optimize model

   model.optimize();

   cout << x.get(GRB_StringAttr_VarName) << " "
   << x.get(GRB_DoubleAttr_X) << endl;
   cout << y.get(GRB_StringAttr_VarName) << " "
   << y.get(GRB_DoubleAttr_X) << endl;
   cout << z.get(GRB_StringAttr_VarName) << " "
   << z.get(GRB_DoubleAttr_X) << endl;

   cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;

   } catch(GRBException e) {
   cout << "Error code = " << e.getErrorCode() << endl;
   cout << e.getMessage() << endl;
   } catch(...) {
   cout << "Exception during optimization" << endl;
   }

}

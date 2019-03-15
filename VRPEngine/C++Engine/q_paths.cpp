#pragma once

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <limits>
#include <iterator>
#include <numeric>
#include "lower_bounds.h"
#include <algorithm>
#include <numeric>
#include <math.h>
#include "prettyprint.hpp"
#include "VRPClass.cpp"

// The struct of results
struct QPathsLim {
   vector<vector<vector<double>>> f;
   vector<vector<vector<double>>> phi;
   vector<vector<vector<int>>> p;
   vector<vector<vector<vector<int>>>> q_route;
   vector<vector<vector<vector<int>>>> q_route_2;
};


// This function returns the q-paths that the q-routes function uses
// to calulate lower bounds
QPathsLim construct_q_paths_lim_(
   VRP& vrp,
   int h,
   vector<vector<double>> &distance_dict,
   string direction,
   int limit
){

   //Calculate lengths
   int len_values = vrp.possible_v[h-vrp.len_N()].values.size();
   int len_N = vrp.len_N();

   //Construct infinity
   double inf = numeric_limits<double>::infinity();
   int inf_int = numeric_limits<int>::max();

   //Initialize the routes
   vector<vector<vector<double>>> f(limit, vector<vector<double>>(len_values, vector<double> (len_N, inf)));
   vector<vector<vector<double>>> phi(limit, vector<vector<double>>(len_values, vector<double> (len_N, inf)));
   vector<vector<vector<int>>> p(limit, vector<vector<int>>(len_values, vector<int> (len_N)));
   vector<vector<vector<vector<int>>>> q_route(limit, vector<vector<vector<int>>>(len_values, vector<vector<int>> (len_N, vector<int>(0))));
   vector<vector<vector<vector<int>>>> q_route_2(limit, vector<vector<vector<int>>>(len_values, vector<vector<int>> (len_N, vector<int>(0))));

   // Initialize some important values
   int truck_capacity = vrp.capacities[h-vrp.len_N()];
   vector<int> &quantities = vrp.quantities;
   vector<int> &values = vrp.possible_v[h-vrp.len_N()].values;
   map<int,int> &values_pos = vrp.possible_v[h-vrp.len_N()].values_pos;
   vector<int> &N = vrp.N;

   // Initialize the routes
   for (int n = 0; n < len_N; n++) {
      int q = vrp.quantities[n];
      if (q <= truck_capacity){
         int l = values_pos[q];
         if (direction == "left"){
            f[0][l][n] = distance_dict[h][N[n]];
         } else {
            f[0][l][n]  = distance_dict[N[n]][h];
         }
         p[0][l][n] = h;
         const int args[] = {h, N[n]};
         q_route[0][l][n].insert(q_route[0][l][n].end(), args, args+2);
      }
   }

   for (int lim = 1; lim < limit; lim++) {
      for(int l = 0; l < len_values; l++) {
         int Q = values[l];
         vector<vector<double>> g(len_N,vector<double>(len_N));
         vector<vector<vector<int>>> g_type(len_N,vector<vector<int>>(len_N,vector<int>(2,0)));
         for (int x_i = 0; x_i < len_N; x_i++){
            int q_p = Q - quantities[x_i];
            if (q_p > 0){
               int l_p = values_pos[q_p];
               for (int x_j = 0; x_j < len_N; x_j++){
                  if (x_i != x_j){
                     if (p[lim-1][l_p][x_j]!=x_i){
                        if (direction == "left"){
                           g[x_i][x_j] = f[lim-1][l_p][x_j] + distance_dict[N[x_j]][N[x_i]];
                        } else {
                           g[x_i][x_j] = f[lim-1][l_p][x_j] + distance_dict[N[x_i]][N[x_j]];
                        }

                        // We save a boolean and the quantity at which g is calc
                        g_type[x_i][x_j][0] = 0;
                        g_type[x_i][x_j][1] = l_p;
                     } else {
                        if (direction == "left"){
                           g[x_i][x_j] = phi[lim-1][l_p][x_j] + distance_dict[N[x_j]][N[x_i]];
                        } else {
                           g[x_i][x_j] = phi[lim-1][l_p][x_j] + distance_dict[N[x_i]][N[x_j]];
                        }
                        g_type[x_i][x_j][0] = 1;
                        g_type[x_i][x_j][1] = l_p;
                     }
                  }
               }
            }
         }

         for (int x_i = 0; x_i < len_N; x_i++){
            int q_p = Q - quantities[x_i];
            if (q_p > 0){
               // Find the minimum
               int arg_min_1 = inf_int;
               int arg_min_2 = inf_int;
               double min_1 = inf;
               double min_2 = inf;
               for (int x_j = 0; x_j < len_N; x_j++){
                  if (x_i!=x_j){
                     double value = g[x_i][x_j];
                     if (value<=min_1){
                        min_2 = min_1;
                        min_1 = value;
                        arg_min_2 = arg_min_1;
                        arg_min_1 = x_j;
                     } else if (value<=min_2) {
                        min_2 = value;
                        arg_min_2 = x_j;
                     }
                  }
               }
               // To fix the bug in case len_N==2
               if (arg_min_2 == inf_int){
                  arg_min_2 = arg_min_1;
                  cout<<"Weird behavior"<<endl;
               }
               p[lim][l][x_i] = arg_min_1;
               f[lim][l][x_i] = min_1;
               phi[lim][l][x_i] = min_2;
               vector<int> &coord = g_type[x_i][arg_min_1];
               vector<int> &coord_2 = g_type[x_i][arg_min_2];

               if (min_1 != inf){
                  q_route[lim][l][x_i] = (coord[0] == 0) ? q_route[lim-1][coord[1]][arg_min_1] : q_route_2[lim-1][coord[1]][arg_min_1];
                  q_route[lim][l][x_i].push_back(x_i);
               }
               if (min_2!= inf){
                  q_route_2[lim][l][x_i] = (coord_2[0] == 0) ? q_route[lim-1][coord_2[1]][arg_min_2] : q_route_2[lim-1][coord_2[1]][arg_min_2];
                  q_route_2[lim][l][x_i].push_back(x_i);
               }



            }
         }
      }
   }

   QPathsLim qpaths;
   qpaths.f = f;
   qpaths.phi = phi;
   qpaths.p = p;
   qpaths.q_route = q_route;
   qpaths.q_route_2 = q_route_2;

   return qpaths;

}

// This function returns the q-paths that the q-routes function uses
// to calulate lower bounds
QPaths construct_q_paths_(
   VRP& vrp,
   int h,
   vector<vector<double>> &distance_dict,
   string direction
){
   //Calculate lengths
   int len_values = vrp.possible_v[h-vrp.len_N()].values.size();
   int len_N = vrp.len_N();


   //Construct infinity
   double inf = numeric_limits<double>::infinity();
   int inf_int = numeric_limits<int>::max();

   //Initialize the routes
   vector<vector<double>> f(len_values, vector<double> (len_N));
   vector<vector<double>> phi(len_values, vector<double> (len_N));
   vector<vector<int>> p(len_values, vector<int> (len_N));
   vector<vector<vector<int>>> q_route(len_values, vector<vector<int>> (len_N, vector<int>(0)));
   vector<vector<vector<int>>> q_route_2(len_values, vector<vector<int>> (len_N, vector<int>(0)));

   // Initialize some important values
   int truck_capacity = vrp.capacities[h-vrp.len_N()];
   vector<int> &quantities = vrp.quantities;
   vector<int> &values = vrp.possible_v[h-vrp.len_N()].values;
   map<int,int> &values_pos = vrp.possible_v[h-vrp.len_N()].values_pos;
   vector<int> &N = vrp.N;


   for(int l = 0; l < len_values; l++) {
      for (int n = 0; n < len_N; n++) {
       f[l][n] = inf;
       phi[l][n] = inf;
      }
   }
   // Initialize the routes
   for (int n = 0; n < len_N; n++) {
      int q = quantities[n];
      if (q <= truck_capacity){
         int l = values_pos[q];
         if (direction == "left"){
            f[l][n] = distance_dict[h][N[n]];
         } else {
         f[l][n]  = distance_dict[N[n]][h];
         }
         p[l][n] = h;
         const int args[] = {h, N[n]};
         q_route[l][n].insert(q_route[l][n].end(), args, args+2);
      }
   }

   for(int l = 0; l < len_values; l++) {

      int Q = values[l];
      vector<vector<double>> g(len_N,vector<double>(len_N));
      vector<vector<vector<int>>> g_type(len_N,vector<vector<int>>(len_N,vector<int>(2,0)));
      for (int x_i = 0; x_i < len_N; x_i++){
         int q_p = Q - quantities[x_i];
         if (q_p > 0){
            int l_p = values_pos[q_p];
            for (int x_j = 0; x_j < len_N; x_j++){
               if (x_i != x_j){
                  if (p[l_p][x_j]!=x_i){
                     if (direction == "left"){
                        g[x_i][x_j] = f[l_p][x_j] + distance_dict[N[x_j]][N[x_i]];
                     } else {
                        g[x_i][x_j] = f[l_p][x_j] + distance_dict[N[x_i]][N[x_j]];
                     }

                     // We save a boolean and the quantity at which g is calc
                     g_type[x_i][x_j][0] = 0;
                     g_type[x_i][x_j][1] = l_p;
                  } else {
                     if (direction == "left"){
                        g[x_i][x_j] = phi[l_p][x_j] + distance_dict[N[x_j]][N[x_i]];
                     } else {
                        g[x_i][x_j] = phi[l_p][x_j] + distance_dict[N[x_i]][N[x_j]];
                     }
                     g_type[x_i][x_j][0] = 1;
                     g_type[x_i][x_j][1] = l_p;
                  }
               }
            }
         }
      }

      for (int x_i = 0; x_i < len_N; x_i++){
         int q_p = Q - quantities[x_i];
         if (q_p > 0){
            // Find the minimum
            int arg_min_1 = inf_int;
            int arg_min_2 = inf_int;
            double min_1 = inf;
            double min_2 = inf;
            for (int x_j = 0; x_j < len_N; x_j++){
               if (x_i!=x_j){
                  double value = g[x_i][x_j];
                  if (value<=min_1){
                     min_2 = min_1;
                     min_1 = value;
                     arg_min_2 = arg_min_1;
                     arg_min_1 = x_j;
                  } else if (value<=min_2) {
                     min_2 = value;
                     arg_min_2 = x_j;
                  }
               }
            }
            // To fix the bug in case len_N==2
            if (arg_min_2 == inf_int){
               arg_min_2 = arg_min_1;
            }
            p[l][x_i] = arg_min_1;
            f[l][x_i] = min_1;
            phi[l][x_i] = min_2;
            vector<int> &coord = g_type[x_i][arg_min_1];
            vector<int> &coord_2 = g_type[x_i][arg_min_2];
            if (min_1 != inf){
               q_route[l][x_i] = (coord[0] == 0) ? q_route[coord[1]][arg_min_1] : q_route_2[coord[1]][arg_min_1];
               q_route[l][x_i].push_back(x_i);
            }
            if (min_2 != inf){
               q_route_2[l][x_i] = (coord_2[0] == 0) ? q_route[coord_2[1]][arg_min_2] : q_route_2[coord_2[1]][arg_min_2];
               q_route_2[l][x_i].push_back(x_i);
            }

         }
      }

   }

   QPaths qpaths;
   qpaths.f = f;
   qpaths.phi = phi;
   qpaths.p = p;
   qpaths.q_route = q_route;
   qpaths.q_route_2 = q_route_2;

   return qpaths;

}


QRoutes construct_q_routes_lim_(
   VRP& vrp,
   int h,
   vector<vector<double>> &distance_dict,
   int limit
){
   //Calculate lengths
   int len_values = vrp.possible_v[h-vrp.len_N()].values.size();
   int len_N = vrp.len_N();

   // Initialize some important values
   int truck_capacity = vrp.capacities[h-vrp.len_N()];
   vector<int> &quantities = vrp.quantities;
   vector<int> &values = vrp.possible_v[h-vrp.len_N()].values;
   map<int,int> &values_pos = vrp.possible_v[h-vrp.len_N()].values_pos;
   vector<int> &N = vrp.N;

   //Construct infinity
   double inf = numeric_limits<double>::infinity();

   QPathsLim qpaths_l = construct_q_paths_lim_(vrp,h,distance_dict,"left",limit);
   QPathsLim qpaths_r = construct_q_paths_lim_(vrp,h,distance_dict,"right",limit);

   vector<vector<vector<double>>>& f_l = qpaths_l.f;
   vector<vector<vector<double>>>& phi_l = qpaths_l.phi;
   vector<vector<vector<int>>>& p_l = qpaths_l.p;
   vector<vector<vector<vector<int>>>>& q_route_l = qpaths_l.q_route;
   vector<vector<vector<vector<int>>>>& q_route_2_l = qpaths_l.q_route_2;

   vector<vector<vector<double>>>& f_r = qpaths_r.f;
   vector<vector<vector<double>>>& phi_r = qpaths_r.phi;
   vector<vector<vector<int>>>& p_r = qpaths_r.p;
   vector<vector<vector<vector<int>>>>& q_route_r = qpaths_r.q_route;
   vector<vector<vector<vector<int>>>>& q_route_2_r= qpaths_r.q_route_2;



   vector<vector<double>> psi(len_values, vector<double> (len_N));
   vector<vector<vector<int>>> psi_route(len_values, vector<vector<int>> (len_N, vector<int>(0)));

   for(int l = 0; l < len_values; l++) {
      for(int n = 0; n < len_N; n++) {
         int min_w = quantities[n];
         // Danger here!!
         int max_w = values[l];
         double min_val = inf;
         vector<int> min_coord(5,-1);
         double val(0);
         int option;
         for (int l_1 = 0; l_1 < len_values; l_1++){
            int q = values[l_1];
            if (q>=min_w && q<=max_w){
               int l_2 = values_pos[values[l] + quantities[n] - q];
               for (int len_1 = 0; len_1 < limit; len_1++){
                  for (int len_2 = 0; len_2 < limit - len_1; len_2++){
                     if (p_l[len_1][l_1][n]!=p_r[len_2][l_2][n] || (p_l[len_1][l_1][n] == p_r[len_2][l_2][n] && p_l[len_1][l_1][n]==h)){
                        val = f_l[len_1][l_1][n]+f_r[len_2][l_2][n];
                        option = 0;
                     } else{
                        if (f_l[len_1][l_1][n]+phi_r[len_2][l_2][n]<phi_l[len_1][l_1][n]+f_r[len_2][l_2][n]){
                           val = f_l[len_1][l_1][n]+phi_r[len_2][l_2][n];
                           option = 1;
                        } else {
                           val = phi_l[len_1][l_1][n]+f_r[len_2][l_2][n];
                           option = 2;
                        }
                     }
                     if (val < min_val){
                        min_val = val;
                        min_coord[0] = option;
                        min_coord[1] = l_1;
                        min_coord[2] = l_2;
                        min_coord[3] = len_1;
                        min_coord[4] = len_2;
                     }
                  }
               }
            }
         }
         psi[l][n] = min_val;
         if (min_coord[0]!=-1){
            if (min_coord[0] == 0){
               vector<int> & A = q_route_l[min_coord[3]][min_coord[1]][n];
               vector<int> & B = q_route_r[min_coord[4]][min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            } else if (min_coord[0] == 1){
               vector<int> & A = q_route_l[min_coord[3]][min_coord[1]][n];
               vector<int> & B = q_route_2_r[min_coord[4]][min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            } else if (min_coord[0] == 2){
               vector<int> & A = q_route_2_l[min_coord[3]][min_coord[1]][n];
               vector<int> & B = q_route_r[min_coord[4]][min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            }
         }
      }
   }

   QRoutes qroutes;
   qroutes.psi = psi;
   qroutes.psi_route = psi_route;

   return qroutes;
}



QRoutes construct_q_routes_(
   VRP& vrp,
   int h,
   vector<vector<double>> &distance_dict
){
   //Calculate lengths
   int len_values = vrp.possible_v[h-vrp.len_N()].values.size();
   int len_N = vrp.len_N();

   // Initialize some important values
   int truck_capacity = vrp.capacities[h-vrp.len_N()];
   vector<int> &quantities = vrp.quantities;
   vector<int> &values = vrp.possible_v[h-vrp.len_N()].values;
   map<int,int> &values_pos = vrp.possible_v[h-vrp.len_N()].values_pos;
   vector<int> &N = vrp.N;

   //Construct infinity
   double inf = numeric_limits<double>::infinity();

   QPaths qpaths_l = construct_q_paths_(vrp,h,distance_dict,"left");
   QPaths qpaths_r = construct_q_paths_(vrp,h,distance_dict,"right");

   vector<vector<double>>& f_l = qpaths_l.f;
   vector<vector<double>>& phi_l = qpaths_l.phi;
   vector<vector<int>>& p_l = qpaths_l.p;
   vector<vector<vector<int>>>& q_route_l = qpaths_l.q_route;
   vector<vector<vector<int>>>& q_route_2_l = qpaths_l.q_route_2;

   vector<vector<double>>& f_r = qpaths_r.f;
   vector<vector<double>>& phi_r = qpaths_r.phi;
   vector<vector<int>>& p_r = qpaths_r.p;
   vector<vector<vector<int>>>& q_route_r = qpaths_r.q_route;
   vector<vector<vector<int>>>& q_route_2_r= qpaths_r.q_route_2;


   vector<vector<double>> psi(len_values, vector<double> (len_N));
   vector<vector<vector<int>>> psi_route(len_values, vector<vector<int>> (len_N, vector<int>(0)));

   for(int l = 0; l < len_values; l++) {
      for(int n = 0; n < len_N; n++) {
         int min_w = quantities[n];
         //some changes here
         int max_w = values[l];
         double min_val = inf;
         vector<int> min_coord(3,-1);
         double val(0);
         int option;
         for (int l_1 = 0; l_1 < len_values; l_1++){
            int q = values[l_1];
            if (q>=min_w && q<=max_w){
               int l_2 = values_pos[values[l] + quantities[n] - q];
               if (p_l[l_1][n]!=p_r[l_2][n] || (p_l[l_1][n] == p_r[l_2][n] && p_l[l_1][n]==h)){
                  val = f_l[l_1][n]+f_r[l_2][n];
                  option = 0;
               } else{
                  if (f_l[l_1][n]+phi_r[l_2][n]<phi_l[l_1][n]+f_r[l_2][n]){
                     val = f_l[l_1][n]+phi_r[l_2][n];
                     option = 1;
                  } else {
                     val = phi_l[l_1][n]+f_r[l_2][n];
                     option = 2;
                  }
               }
               if (val < min_val){
                  min_val = val;
                  min_coord[0] = option;
                  min_coord[1] = l_1;
                  min_coord[2] = l_2;
               }
            }
         }
         psi[l][n] = min_val;
         if (min_coord[0]!=-1){
            if (min_coord[0] == 0){
               vector<int> & A = q_route_l[min_coord[1]][n];
               vector<int> & B = q_route_r[min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            } else if (min_coord[0] == 1){
               vector<int> & A = q_route_l[min_coord[1]][n];
               vector<int> & B = q_route_2_r[min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            } else if (min_coord[0] == 2){
               vector<int> & A = q_route_2_l[min_coord[1]][n];
               vector<int> & B = q_route_r[min_coord[2]][n];
               psi_route[l][n].reserve( A.size() + B.size() - 1 ); // preallocate memory
               psi_route[l][n].insert( psi_route[l][n].end(), A.begin(), A.end() );
               psi_route[l][n].insert( psi_route[l][n].end(), B.rbegin() + 1, B.rend() ); // insert backwards
            }
         }
      }
   }

   QRoutes qroutes;
   qroutes.psi = psi;
   qroutes.psi_route = psi_route;

   return qroutes;
}


LB_GENPATH path_lower_bound(
   int h,
   VRP& vrp,
   vector<vector<double>> &distance_dict,
   string direction,
   int limit
){

   double inf = numeric_limits<double>::infinity();
   int inf_int = numeric_limits<int>::max();




   int truck_capacity = vrp.capacities[h-vrp.len_N()];
   vector<int>& N = vrp.N;
   vector<int> &values = vrp.possible_v[h-vrp.len_N()].values;
   map<int,int> &values_pos = vrp.possible_v[h-vrp.len_N()].values_pos;
   vector<int> &quantities = vrp.quantities;

   int len_values = vrp.possible_v[h-vrp.len_N()].values.size();
   int len_N = vrp.len_N();



   QPaths qpaths;
   if (limit != inf_int){
      // Generate the qpaths_lim
      QPathsLim qpaths_lim = construct_q_paths_lim_(vrp,
        h,
        distance_dict,
        direction,
        limit
      );
      // Take them and transform them into qpaths
      //Initialize the routes
      vector<vector<double>> f(len_values, vector<double> (len_N));
      vector<vector<double>> phi(len_values, vector<double> (len_N));
      vector<vector<int>> p(len_values, vector<int> (len_N));
      vector<vector<vector<int>>> q_route(len_values, vector<vector<int>> (len_N, vector<int>(0)));
      vector<vector<vector<int>>> q_route_2(len_values, vector<vector<int>> (len_N, vector<int>(0)));

      for (int l = 0; l < len_values; l++) {
         for (int end = 0; end<len_N; end++){
            int argmin = -1;
            double min_val = inf;
            // Find minimum
            for (int lim = 0; lim < limit; lim++) {
               if (qpaths_lim.f[lim][l][end] <= min_val){
                  argmin = lim;
                  min_val = qpaths_lim.f[lim][l][end];
               }
            }
            f[l][end] = min_val;

            p[l][end] = qpaths_lim.p[argmin][l][end];
            q_route[l][end] = qpaths_lim.q_route[argmin][l][end];
            // Find the phi
            int argmin2 = -1;
            double min_val2 = inf;
            //whether the found route is an f or a phi
            int option = -1;
            for (int lim = 0; lim < limit; lim++) {
               if (qpaths_lim.p[lim][l][end] == p[l][end]){
                  if (qpaths_lim.phi[lim][l][end] <= min_val2){
                     argmin2 = lim;
                     min_val2 = qpaths_lim.phi[lim][l][end];
                     option = 0;
                  }
               } else {
                  if (qpaths_lim.f[lim][l][end] <= min_val2){
                     argmin2 = lim;
                     min_val2 = qpaths_lim.f[lim][l][end];
                     option = 1;
                  }
               }
            }
            phi[l][end] = min_val2;
            if (option == 0){
               q_route_2[l][end] = qpaths_lim.q_route_2[argmin2][l][end];
            } else if (option == 1){
               q_route_2[l][end] = qpaths_lim.q_route[argmin2][l][end];
            }
         }
      }
      qpaths.f = f;
      qpaths.phi = phi;
      qpaths.p = p;
      qpaths.q_route = q_route;
      qpaths.q_route_2 = q_route_2;

  } else{
      qpaths = construct_q_paths_(
         vrp,
         h,
         distance_dict,
         direction
      );
   }



   //Debug
   // if (direction == "right"){
   //    int q = 90;
   //    int p = 62;
   //    cout<<"qpath"<<endl;
   //    cout<<qpaths.q_route[values_pos[q]][p]<<endl;
   //    cout<<qpaths.f[values_pos[q]][p]<<endl;
   //    double real_cost = 0;
   //    for (int s = 0; s < (int) qpaths.q_route[values_pos[q]][p].size() - 1; s++){
   //       if (direction=="right"){
   //          real_cost  += distance_dict[qpaths.q_route[values_pos[q]][p][s+1]][qpaths.q_route[values_pos[q]][p][s]];
   //       } else {
   //          real_cost  += distance_dict[qpaths.q_route[values_pos[q]][p][s]][qpaths.q_route[values_pos[q]][p][s+1]];
   //       }
   //    }
   //    cout<<real_cost<<endl;
   // }

   // Define infinity
   double infinity = numeric_limits<double>::infinity();


   // Initialize values
   vector<vector<double>> F(len_values,vector<double>(len_N,infinity));
   vector<vector<double>> G(len_values,vector<double>(len_N,infinity));
   vector<vector<int>> X(len_values,vector<int>(len_N,-1));
   vector<vector<vector<int>>> min_q_path(len_values,vector<vector<int>>(len_N));
   vector<vector<vector<int>>> min_q_path_2(len_values,vector<vector<int>>(len_N));


   // The F function takes as argument the farmer
   for (int i = 0; i < len_N; i++){
      int q_i = quantities[i];
      int q_lb = values_pos[q_i];
      int q_ub = values_pos[truck_capacity];
      for (int pos_q = q_lb; pos_q <= q_ub; pos_q++){
         int qp_lb = q_lb;
         int qp_ub = values_pos[truck_capacity - values[pos_q] + q_i];
         double min_F = infinity;
         int arg_min_F = -1;
         for (int qp = qp_lb; qp <= qp_ub; qp++){
            if ((qpaths.f)[qp][i] <= min_F){
               min_F = (qpaths.f)[qp][i];
               arg_min_F = qp;
            }
         }

         F[pos_q][i] = min_F;
         F[pos_q][i] = (qpaths.f)[arg_min_F][i];
         X[pos_q][i] = (qpaths.p)[arg_min_F][i];
         min_q_path[pos_q][i] = (qpaths.q_route)[arg_min_F][i];
         double min_G = infinity;
         int arg_min_G = -1;
         bool used_f = true;
         //cout<<qp_lb<<","<<qp_ub<<endl;
         for (int qp = qp_lb; qp <= qp_ub; qp++){
            if (qpaths.p[qp][i] == X[pos_q][i]){
               if ((qpaths.phi)[qp][i] <= min_G){
                  min_G = (qpaths.phi)[qp][i];
                  arg_min_G = qp;
                  used_f = false;
               }
            } else {
               if ((qpaths.f)[qp][i] <= min_G){
                  min_G = (qpaths.f)[qp][i];
                  arg_min_G = qp;
                  used_f = true;
               }
            }
         }
         G[pos_q][i] = min_G;
         if (used_f){
            min_q_path_2[pos_q][i] = (qpaths.q_route)[arg_min_G][i];
         } else {
            min_q_path_2[pos_q][i] = (qpaths.q_route_2)[arg_min_G][i];
         }

      }

   }

   LB_GENPATH functions;
   functions.F = F;
   functions.X = X;
   functions.G = G;
   functions.min_q_path = min_q_path;
   functions.min_q_path_2 = min_q_path_2;

   //Debug
   // if (direction == "right"){
   //    cout<<"min_q_path"<<endl;
   //    cout<<min_q_path[values_pos[10]][62]<<endl;
   //    cout<<F[values_pos[10]][62]<<endl;
   // }


   return functions;

}

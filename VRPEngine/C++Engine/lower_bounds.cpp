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
#include "q_paths.cpp"



template class std::vector<double>;
template class std::vector<std::vector<double>>;




using namespace std;

// This function takes a disance matrix and a vector on vertex penalties on the farmers, and then it penalizes
vector<vector<double>> penalized_matrix(
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


vector<vector<vector<double>>> calc_reduced_truck_distances(
   vector<vector<vector<double>>> truck_distances,
   vector<double> lamb,
   vector<double> mu
){
   vector<vector<vector<double>>> reduced_truck_distances;
   int len_H = truck_distances.size();
   for(int i = 0; i<len_H; i++){
      reduced_truck_distances.push_back(reduced_cost_matrix(truck_distances[i], lamb, mu));
   }
   return reduced_truck_distances;
}

vector<vector<double>> reduced_cost_matrix(
   vector<vector<double>> truck_distance,
   vector<double> lamb,
   vector<double> mu
){
   vector<double> penalties;
   penalties.reserve( lamb.size() + mu.size() ); // preallocate memory
   penalties.insert( penalties.end(), lamb.begin(), lamb.end() );
   penalties.insert( penalties.end(), mu.begin(), mu.end() );
   for (int i = 0; i < (int) truck_distance.size(); i++){
      for (int j = 0; j < (int) truck_distance.size(); j++){
         truck_distance[i][j] -= (penalties[i]/2.0 + penalties[j]/2.0);
      }
   }
   return truck_distance;
}


void DualSolution::sol_calc_reduced_truck_distances(vector<vector<vector<double>>> truck_distances){
   reduced_truck_distances = calc_reduced_truck_distances(truck_distances,u,v);
}



void DualSolution::initialize_routes(int len_H){
   routes = vector<list<SimpleRoute>>(len_H);
}




// Write lower bounds
LowerBound lower_bound_(
   VRP &vrp,
   vector<vector<vector<double>>> &red_truck_distances,
   vector<double> mu,
   vector<double> lamb,
   int limit
){

   vector<int>& H = vrp.H;
   vector<int>& capacities = vrp.capacities;
   vector<int>& N = vrp.N;
   vector<int>& quantities = vrp.quantities;
   vector<int>& n_trucks = vrp.n_trucks;

   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();


   //The vectors of minima
   vector<vector<double>> b(len_N, vector<double>(len_H));
   vector<vector<int>> val(len_N, vector<int>(len_H));
   vector<vector<vector<int>>> b_routes(len_N, vector<vector<int>>(len_H, vector<int>(0)));
   for (int h = 0; h < len_H; h++){
      int truck_capacity = capacities[h];

      PossibleValues possible = possible_values(quantities, truck_capacity);

      QRoutes qroutes;
      bool old_code = true;
      if (limit == numeric_limits<int>::max()){
         if (old_code){
            qroutes = construct_q_routes_(H[h], truck_capacity, N, red_truck_distances[h], possible.values, possible.values_pos, quantities);
         } else {
            qroutes = construct_q_routes_(vrp, H[h], red_truck_distances[h]);
         }
      } else {
         if (!old_code){
            qroutes = construct_q_routes_lim_(vrp, H[h], red_truck_distances[h], limit);
         } else {
            qroutes = construct_q_routes_(H[h], truck_capacity, N, red_truck_distances[h], possible.values, possible.values_pos, quantities);
         }
      }


      // We find the minimum l

      int len_values = possible.values.size();
      for (int n = 0; n < len_N; n++){
         vector<double> b_values(len_values);

         for (int l = 0; l < len_values; l++){
            //b_values[l] = (qroutes.psi[l][n] - mu[h]) * (double)quantities[n]/(double)possible.values[l];
            b_values[l] = (qroutes.psi[l][n]) * (double)quantities[n]/(double)possible.values[l];
         }

         std::vector<double>::iterator l_it = std::min_element(b_values.begin(), b_values.end());
         int l_min = std::distance(b_values.begin(), l_it);
         b[n][h] = b_values[l_min];
         b_routes[n][h] = qroutes.psi_route[l_min][n];
         val[n][h] = possible.values[l_min];


      }


   }


   // The second vector of minima
   vector<double> b_min(len_N);
   vector<int> val_min(len_N);
   vector<int> h_min(len_N);
   vector<vector<int>> b_min_routes(len_N, vector<int>(0));
   for (int n = 0; n < len_N; n++){
      std::vector<double>::iterator h_it = std::min_element(b[n].begin(), b[n].end());
      int h_m = std::distance(b[n].begin(), h_it);
      b_min[n] = b[n][h_m];
      b_min_routes[n] = b_routes[n][h_m];
      val_min[n] = val[n][h_m];
      h_min[n] = h_m;
   }

   //Calculate number of visits to each node
   vector<vector<int>> visits(len_N,vector<int>(len_N,0));
   for (int n = 0; n < len_N; n++){
      for(vector<int>::iterator it_route = b_min_routes[n].begin();  it_route!=b_min_routes[n].end(); ++it_route){
         if (*it_route < len_N){
            visits[n][*it_route] += 1;
         }
      }
   }


   // Construct theta
   vector<double> theta(len_N,0);
   for (int j = 0; j<len_N; j++){
      for (int i=0; i<len_N; i++){
         theta[j] += ((double) quantities[i])/((double) val_min[i])*(double)visits[i][j];
      }
   }
   // Construct rho
   vector<double> rho(len_H,0);
   for (int n=0; n<len_N; n++){
      rho[h_min[n]] += ((double) quantities[n])/((double) val_min[n]);
   }

   // Construct dual variables
   vector<double> u(len_N,0);
   for (int n=0; n<len_N; n++){
      u[n] = b_min[n] + lamb[n];
   }

   double z_lb = 0;
   for (int n=0; n<len_N; n++){
      z_lb += u[n];
   }
   for (int h=0; h<len_H; h++){
      z_lb += mu[h]*n_trucks[h];
   }


   LowerBound lb;
   lb.z_lb = z_lb;
   lb.u = u;
   lb.theta = theta;
   lb.rho = rho;

   return lb;

}


// This function returns the q-paths that the q-routes function uses
// to calulate lower bounds
QPaths construct_q_paths_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> &distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities,
   string direction
){
   //Calculate lengths
   int len_values = values.size();
   int len_N = N.size();

   //Construct infinity
   double inf = numeric_limits<double>::infinity();
   int inf_int = numeric_limits<int>::max();

   //Initialize the routes
   vector<vector<double>> f(len_values, vector<double> (len_N));
   vector<vector<double>> phi(len_values, vector<double> (len_N));
   vector<vector<int>> p(len_values, vector<int> (len_N));
   vector<vector<vector<int>>> q_route(len_values, vector<vector<int>> (len_N, vector<int>(0)));
   vector<vector<vector<int>>> q_route_2(len_values, vector<vector<int>> (len_N, vector<int>(0)));

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
         // We need the quantity len_N to be >1 to avoid bugs
         if (q_p > 0 && len_N > 1){
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
            q_route[l][x_i] = (coord[0] == 0) ? q_route[coord[1]][arg_min_1] : q_route_2[coord[1]][arg_min_1];
            q_route_2[l][x_i] = (coord_2[0] == 0) ? q_route[coord_2[1]][arg_min_2] : q_route_2[coord_2[1]][arg_min_2];
            q_route[l][x_i].push_back(x_i);
            q_route_2[l][x_i].push_back(x_i);

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

// Constructs functions to calculate lower bounds of a path, that are used later
// by the GENROUTE function
LB_GENPATH path_lower_bound(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> &distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities,
   string direction
){
   //cout<<direction<<endl;
   //Construct the qpaths
   QPaths qpaths = construct_q_paths_(h, truck_capacity, N, distance_dict, values, values_pos, quantities, direction);
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
   int len_N = N.size();
   int len_values = values.size();

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


QRoutes construct_q_routes_(
   int h,
   int truck_capacity,
   vector<int> N,
   vector<vector<double>> &distance_dict,
   vector<int> values,
   map<int,int> values_pos,
   vector<int> quantities
){


   QPaths qpaths_l = construct_q_paths_(h,truck_capacity,N,distance_dict,values,values_pos,quantities,"left");
   QPaths qpaths_r = construct_q_paths_(h,truck_capacity,N,distance_dict,values,values_pos,quantities,"right");

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




   //Calculate lengths
   int len_values = values.size();
   int len_N = N.size();

   //Construct infinity
   double inf = numeric_limits<double>::infinity();

   vector<vector<double>> psi(len_values, vector<double> (len_N));
   vector<vector<vector<int>>> psi_route(len_values, vector<vector<int>> (len_N, vector<int>(0)));

   for(int l = 0; l < len_values; l++) {
      for(int n = 0; n < len_N; n++) {
         int min_w = quantities[n];
         int max_w = values[l] + quantities[n];
         double min_val = inf;
         vector<int> min_coord(3,-1);
         double val(0);
         int option;
         for (int l_1 = 0; l_1 < len_values; l_1++){
            int q = values[l_1];
            if (q>=min_w && q<max_w){
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

// Given a set of routes, calculate the best lower bounds that they generate
// We initialize at mu and lamb
DualSolution lower_bound_optimizer_M1(
   int iterations,
   double z_ub,
   double epsilon,
   VRP &vrp,
   int limit
){

   vector<int>& H = vrp.H;
   vector<int>& capacities = vrp.capacities;
   vector<int>& N = vrp.N;
   vector<int>& quantities = vrp.quantities;
   vector<int>& n_trucks = vrp.n_trucks;

   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();
   // Define infinity
   double infinity = numeric_limits<double>::infinity();
   // Define u
   vector<double> u(len_N);

   // Initialize mu and lambda
   vector<double> lamb(len_N, 0);
   vector<double> mu(len_H, 0);

   // Vectors to store the optimal values
   vector<double> lamb_opt(len_N);
   vector<double> u_opt(len_N);
   vector<double> v_opt(len_H);
   // Here we store the values of the iterations
   vector<double> values;
   double max_val = -infinity;

   for (int iteration = 0; iteration<iterations; iteration++){


      vector<vector<vector<double>>> red_truck_distances = calc_reduced_truck_distances(vrp.truck_distances, lamb, mu);

      // We pass these routes to the algorithm that calculates the lower bound
      LowerBound lb = lower_bound_(vrp, red_truck_distances, mu, lamb, limit);
      cout<<lb.z_lb<<endl;
      //cout<<lamb<<endl;
      //cout<<mu<<endl;




      // Check if the lower bound that we get improves
      if (lb.z_lb > max_val){
         max_val = lb.z_lb;
         u_opt = lb.u;
         v_opt = mu;
         lamb_opt = lamb;
      }
      values.push_back(lb.z_lb);

      // We calculate g for the step of the algorithm
      double g_den_1 = 0;
      double g_den_2 = 0;
      for (int i = 0; i< len_N; i++){
         g_den_1 += pow(lb.theta[i]-1.0,2);
      }
      for (int i = 0; i< len_H; i++){
         g_den_2 += pow(lb.rho[i]-n_trucks[i],2);
      }
      double g = (z_ub - lb.z_lb)/(g_den_1 + g_den_2);
      // Update lambda and mu
      for (int i = 0; i<len_N; i++){
         lamb[i] = lamb[i] - epsilon*g*(lb.theta[i]-1.0);
      }
      for (int i = 0; i<len_H; i++){
         mu[i] = min(mu[i] - epsilon*g*(lb.rho[i]-n_trucks[i]),0.0);
      }

      // Check that we are not getting exploting reduced variables
      double explotion = 0;
      for (int i = 0; i< len_N; i++){
         explotion+=fabs(lamb[i]);
      }
      if (explotion > pow(10,16)){
         cout<<"Diverging reduced variables"<<endl;
         break;
      }

      // We update epsilon
      if ((int) values.size() >= 7){
         // Calculate changes
         vector<double> grad;
         for (int i = 0; i < (int)values.size()-1; i++){
            grad.push_back(values[i+1] - values[i]);
         }
         // Calculate jumps
         vector<int> jumps;
         for (int i = 0; i < (int)grad.size()-1; i++){
            jumps.push_back(signbit(grad[i+1])!=signbit(grad[i]));
         }
         // If too many jumps reduce epsilon
         int n_jumps = 0;
         for (int i = (int)jumps.size()-5; i<(int)jumps.size(); i++ ){
            n_jumps += jumps[i];
         }
         if (n_jumps >= 3){
            epsilon = epsilon/1.5;
            cout<<"New epsilon "<<epsilon<<endl;
            values.clear();
         }
      }
      // Check if we have reached a zero gradient
      bool zero_gradient = true;
      for (int i = 0; i < len_N; i++){
         if (lb.theta[i] != 1.0){
            zero_gradient = false;
            break;
         }
      }
      for (int i = 0; i < len_H; i++){
         if (lb.rho[i] != n_trucks[i]){
            zero_gradient = false;
            break;
         }
      }
      if (zero_gradient){
         cout<<"Reached zero gradient"<<endl;
         break;
      }
   }

   DualSolution new_bound;
   new_bound.z_lb = max_val;
   new_bound.u = u_opt;
   new_bound.v = v_opt;
   new_bound.lamb = lamb_opt;


   return new_bound;
}

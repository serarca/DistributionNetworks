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
#include <tuple>
#include <list>
#include "lower_bounds.h"
#include "baldacci.h"

#include <cmath>
#include "prettyprint.hpp"
#include <bitset>
#include "VRPClass.cpp"
#include "utils/save_dual_solution.cpp"

#include <chrono>



void print_set(set<tuple<int,int>> s){
   for (auto n:s)
      cout<<"("<<get<0>(n)<<","<<get<1>(n)<<")";
   cout<<endl;
}

void print_path(Path p){
   cout<<"\t Path: ";
   for (auto n:p.path){
      cout<<n<<" ";
   }
   cout<<"Set: ";
   for (int i = 0; i < (int)p.nodes.size(); i++){
      if (p.nodes[i] == 1)
         cout<<i<<" ";
   }
   cout<<"Cost: "<<p.cost<<" ";
   cout<<"Load: "<<p.load<<" ";
   cout<<"Length: "<<p.length<<" ";
   cout<<"LB Path: "<<p.lb_path<<" ";
   cout<<"LB Load: "<<p.lb_load<<" ";
   cout<<"Real Cost: "<<p.real_cost<<" ";
   cout<<"Real LB: "<<p.real_lb<<" ";

   cout<<"Lower Bound: "<<p.lower_bound<<" "<<endl;
}

void print_route(Route r){
   cout<<"\t Route_l: ";
   for (auto n:(r.path_l->path)){
      cout<<n<<" ";
   }
   cout<<"\t Route_r: ";
   for (auto n:(r.path_r->path)){
      cout<<n<<" ";
   }
   cout<<"Set: ";
   for (int i = 0; i < (int)r.nodes.size(); i++){
      if (r.nodes[i] == 1)
         cout<<i<<" ";
   }
   cout<<"Cost: "<<r.cost<<" ";
   cout<<"Load: "<<r.load<<" ";
   cout<<"Median: "<<r.median<<" "<<endl;
}

void print_sroute(SimpleRoute r){
   cout<<"\t Route: ";
   for (auto n:(r.path)){
      cout<<n<<" ";
   }
   cout<<"Cost: "<<r.cost<<" ";
   cout<<"Geo Cost: "<<r.geo_cost<<" ";
   cout<<"Load: "<<r.load<<" ";
   cout<<"Median: "<<r.median<<" ";
   cout<<"Lower_bounds: "<<r.l_lb<<", "<<r.r_lb<<endl;
}

void print_Paths(vector<list<Path>> paths){
   int i = 0;
   for (auto end:paths){
      cout<<"End Node: "<<i;
      cout<<"Length: "<<end.size()<<endl;
      i+=1;

      for (auto p:end){
         print_path(p);
      }
   }
}

void print_paths(list<Path> end){

      cout<<"Length: "<<end.size()<<endl;

      for (auto p:end){
         print_path(p);
      }
}

void print_Routes(vector<list<Route>> routes){
   int i = 0;
   for (auto end:routes){
      cout<<"End Node: "<<i;
      cout<<"Length: "<<end.size()<<endl;
      i+=1;

      for (auto r:end){
         print_route(r);
      }
   }
}

void print_sRoutes(list<SimpleRoute> end){
      cout<<"Length: "<<end.size()<<endl;

      for (auto r:end){
         print_sroute(r);
      }
}

void print_SRoutes(vector<list<SimpleRoute>> routes){
   int i = 0;
   for (auto end:routes){
      cout<<"End Node: "<<i;
      cout<<"Length: "<<end.size()<<endl;
      i+=1;

      for (auto r:end){
         print_sroute(r);
      }
   }
}

void p_v(vector<double> vec){
   cout<<vec<<endl;
}

void p_v_v(vector<vector<double>> vec){
   for (int i = 0; i< (int)vec.size(); i++){
      cout<<i<<":"<<vec[i]<<endl;
   }

}

void p_v_v_v(vector<vector<vector<int>>> vec){
   for (int i = 0; i< (int)vec.size(); i++){
      for (int j = 0; j< (int)vec[i].size(); j++)
         cout<<i<<":"<<j<<":"<<vec[i][j]<<endl;
   }

}




using GenPath = vector<list<Path>>;
using GenRoute = vector<list<SimpleRoute>>;

bool compare_routes (Route i, Route j) { return (i.cost<j.cost); }

vector<list<Path>> GENPATH(
   int Delta,
   double gamma,
   int h,
   int capacity,
   vector<int> N,
   vector<int> quantities,
   vector<vector<double>> &distance_dict,
   string direction,
   bool &terminated,
   double &gamma_guarantee,
   double max_length
){

   auto start5 = std::chrono::high_resolution_clock::now();

   int len_N = N.size();

   PossibleValues pv= possible_values(quantities, capacity);

   // Generate lower bound paths of the opposite direction
   cout<<"Generating lower bounds"<<endl;
   LB_GENPATH Path_lb;
   if (direction == "left")
   {
      Path_lb = path_lower_bound(h,capacity,N,distance_dict,pv.values,pv.values_pos,quantities,"right");
   } else {
      Path_lb = path_lower_bound(h,capacity,N,distance_dict,pv.values,pv.values_pos,quantities,"left");
   }
   cout<<"Finished lower bounds"<<endl;

   vector<list<Path>> P(N.size() + 1, list<Path>(0));
   //vector<list<Path>> T(N.size() + 1, list<Path>(0));

   vector<vector<list<Path>>> P_q(capacity + 1,vector<list<Path>>(N.size() + 1, list<Path>(0)));
   vector<vector<list<Path>>> T_q(capacity + 1,vector<list<Path>>(N.size() + 1, list<Path>(0)));

   Path init;
   init.path.push_front(h);
   init.cost = 0;
   init.lower_bound = 0;
   init.load = 0;
   init.end = h;
   init.length = 0;
   T_q[0][N.size()].push_front(init);

   // Define infinity
   double inf = numeric_limits<double>::infinity();
   int count_paths = 0;


   int cum1 = 0;
   int cum2 = 0;
   int cum3 = 0;
   int cum4 = 0;


   while (true){

      auto start1 = std::chrono::high_resolution_clock::now();
      //
      // if (count_paths % 100 ==0 )
      // {
      //    cout<< count_paths<<endl;
      // }
      // Check the route with the smallest cost
      int arg_min = -1;
      double cost_min = inf;
      int qu_min = -1;
      // No typo here with the <=
      for(int qu = 0; qu <= capacity; qu++)
      {
         for(int i = 0; i <= len_N; i++){
            if (T_q[qu][i].size() > 0){
               double cost = T_q[qu][i].front().lower_bound;
               if (cost<cost_min){
                  cost_min = cost;
                  arg_min = i;
                  qu_min = qu;
               }
            }
         }
      }
      // Break if no more routes
      if (arg_min == -1){
         if (max_length == numeric_limits<double>::infinity()){
            gamma_guarantee = gamma;
            cout<<"     Paths Reached Gamma: "<<gamma<<endl;
            // int total_paths = 0;
            // for (auto p:P){
            //    total_paths+= p.size();
            // }
            // cout<<"        Total No. of Paths: "<<total_paths<<endl;
         }
         // for (int s = 0; s < N.size() + 1; s++){
         //    cout<<"           Farmer "<<s<<":"<<P[s].size()<<endl;
         // }
         // for (auto i_1:P[6]){
         //    print_path(i_1);
         // }
         break;
      }


      // Check the first route and add it to P
      P_q[qu_min][arg_min].splice(P_q[qu_min][arg_min].end(),T_q[qu_min][arg_min],T_q[qu_min][arg_min].begin());
      count_paths += 1;

      // Break if too many paths
      if (count_paths==Delta){
         if (max_length == numeric_limits<double>::infinity()){
            terminated = false;
            gamma_guarantee = cost_min;
            cout<<"     Paths Reached Delta, best LB: "<<cost_min<<endl;
            // int total_paths = 0;
            // for (auto p:P){
            //    total_paths+= p.size();
            // }
            // cout<<"        Total No. of Paths: "<<total_paths<<endl;
         }
         // for (int s = 0; s < N.size() + 1; s++){
         //    cout<<"           Farmer "<<s<<":"<<P[s].size()<<endl;
         // }
         // for (auto i_1:P[6]){
         //    print_path(i_1);
         // }
         break;
      }

      // Extract the element
      Path p_star = P_q[qu_min][arg_min].back();

      //Some debugging
      /*
      Path suspicious;
      suspicious.path.push_back(2);
      suspicious.path.push_back(1);
      if (p_star.path == suspicious.path && direction == "left"){
         cout<<"found it again"<<endl;
      }
      */

      // If path too long, go to the next one

      if ((double) p_star.length >= max_length)
         continue;

      // If path violates capacity, go to the next one
      if ((double) p_star.load > ((double) capacity)/2.0)
         continue;
      // Add new paths
      cum1 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start1).count();


      for (int i = 0; i < len_N; i++){


         // Check node is not in path
         if (p_star.nodes[i] == 0){
            auto start2 = std::chrono::high_resolution_clock::now();

            Path new_path;
            new_path.path = p_star.path;
            new_path.path.push_back(i);
            new_path.nodes = p_star.nodes;
            new_path.nodes[i] = 1;
            new_path.load = p_star.load + quantities[i];
            new_path.end = N[i];
            new_path.length = p_star.length + 1;
            if (direction == "left"){
               new_path.cost = p_star.cost + distance_dict[p_star.end][i];
            }
            if (direction == "right"){
               new_path.cost = p_star.cost + distance_dict[i][p_star.end];
            }

            // We will do some debugging here
            /*
            if (new_path.path == suspicious.path && direction == "left"){
               cout<<"found it"<<endl;
               cout<< Path_lb.X[pv.values_pos[new_path.load]][i] << endl;
               cout<< Path_lb.F[pv.values_pos[new_path.load]][i] << endl;
               cout<< Path_lb.G[pv.values_pos[new_path.load]][i] << endl;
               cout<<Path_lb.min_q_path[pv.values_pos[new_path.load]][i]<<endl;
               cout<<Path_lb.min_q_path_2[pv.values_pos[new_path.load]][i]<<endl;
               p_v_v(distance_dict);
            }
            */



            // Calculate a lower bound
            // Check the next node in the remaining q-path
            int previous_node = Path_lb.X[pv.values_pos[new_path.load]][i];
            bool in_path = false;
            for (auto node:new_path.path){
               if((node == previous_node) && (node<len_N)){
                  in_path = true;
                  break;
               }
            }

            double remaining_cost;

            //vector<int> lb_path;

            if (!in_path){
               remaining_cost = Path_lb.F[pv.values_pos[new_path.load]][i];
               //lb_path = Path_lb.min_q_path[pv.values_pos[new_path.load]][i];
            } else {
               remaining_cost = Path_lb.G[pv.values_pos[new_path.load]][i];
               //lb_path = Path_lb.min_q_path_2[pv.values_pos[new_path.load]][i];
            }
            new_path.lower_bound = new_path.cost + remaining_cost;
            //new_path.lb_path = lb_path;
            ////Calculate load of path
            //new_path.lb_load = 0;

            // for (int s = 1; s<(int) new_path.lb_path.size() - 1; s++){
            //    int n = new_path.lb_path[s];
            //    new_path.lb_load += quantities[n];
            // }
            //
            // //Calculate real lower bound to debug
            // new_path.real_cost = 0;
            // new_path.real_lb = 0;
            // if(new_path.load == 10 && new_path.end == 62)
            // {
            //    cout<<direction<<endl;
            // }
            // for (auto it = new_path.path.begin(); std::next(it,1)!=new_path.path.end(); ++it){
            //    int n = *it;
            //    int n_next = *std::next(it,1);
            //    if (direction=="left"){
            //       new_path.real_cost += distance_dict[n][n_next];
            //       if(new_path.load == 10 && new_path.end == 62)
            //       {
            //          cout<<n<<','<<n_next<<": "<<distance_dict[n][n_next]<<endl;
            //
            //       }
            //    } else {
            //       new_path.real_cost += distance_dict[n_next][n];
            //       if(new_path.load == 10 && new_path.end == 62){
            //          cout<<n<<','<<n_next<<": "<<distance_dict[n_next][n]<<endl;
            //       }
            //    }
            // }
            // new_path.real_lb = new_path.real_cost;
            // for (int s = 0; s < (int) new_path.lb_path.size() - 1; s++){
            //    if (direction=="left"){
            //       new_path.real_lb  += distance_dict[new_path.lb_path[s+1]][new_path.lb_path[s]];
            //       if(new_path.load == 10 && new_path.end == 62){
            //          cout<<distance_dict[new_path.lb_path[s+1]][new_path.lb_path[s]]<<endl;
            //       }
            //    } else {
            //       new_path.real_lb  += distance_dict[new_path.lb_path[s]][new_path.lb_path[s+1]];
            //       if(new_path.load == 10 && new_path.end == 62)
            //       {
            //          cout<<distance_dict[new_path.lb_path[s]][new_path.lb_path[s+1]]<<endl;
            //       }
            //    }
            // }
            //


            /*
            //Debug here
            if (new_path.path == suspicious.path && direction == "left"){
               cout<<new_path.lower_bound<<endl;
               cout<< new_path.load << endl;
               cout<< new_path.cost << endl;
            }
            */
            //Check if the new path has a cost too high
            if (new_path.lower_bound > gamma)
               continue;
            //Check if the new path has a load too high
            if (new_path.load > capacity)
               continue;

            cum2 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start2).count();
            auto start3 = std::chrono::high_resolution_clock::now();

            //Check if this new path is dominated by any path in P
            bool dominated = false;
            for (auto p = P_q[new_path.load][i].begin(); p!=P_q[new_path.load][i].end(); ++p){
                if ((p->end == new_path.end) && (p->cost <= new_path.cost) && (p->nodes == new_path.nodes)){
                   dominated = true;
                   break;
                }
            }
            cum3 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start3).count();
            auto start4 = std::chrono::high_resolution_clock::now();

            /*
            //Debug
            if (new_path.path == suspicious.path && direction == "left"){
               cout<<dominated<<endl;
            }
            */
            if (dominated)
                continue;

            // We check that the path is not dominated, and erase those paths
            // which are dominated
            auto p = T_q[new_path.load][i].begin();
            list<Path>::iterator p_insert;
            bool found_insertion = false;
            while(p!=T_q[new_path.load][i].end()){
               if ((p->end == new_path.end) && (p->cost <= new_path.cost) && (p->nodes == new_path.nodes)){
                  dominated = true;
                  break;
               }
               if ((p->end == new_path.end) && (p->cost > new_path.cost) && (p->nodes == new_path.nodes)){
                  p=T_q[new_path.load][i].erase(p);
                  continue;
               }
               if ((!found_insertion) && (p->lower_bound > new_path.lower_bound)){
                  found_insertion = true;
                  p_insert = p;
               }
               ++p;
            }
            if (!found_insertion){
               p_insert = T_q[new_path.load][i].end();
            }
            T_q[new_path.load][i].insert(p_insert,new_path);
            cum4 += std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start4).count();
         }
      }
   }

   /*
   // We add the one node paths if they are not added already, and we add them
   // at the beginning
   for (int i = 0; i < len_N; i++){
      bool add = true;
      Path new_path;
      new_path.path.push_back(h);
      new_path.path.push_back(i);
      new_path.nodes.insert(i);
      new_path.load = quantities[i];
      new_path.end = N[i];
      if (direction == "left"){
         new_path.cost = distance_dict[h][i];
         new_path.lower_bound = new_path.cost;
      }
      if (direction == "right"){
         new_path.cost = distance_dict[i][h];
         new_path.lower_bound = new_path.cost;
      }
      for (auto p:P[i]){
         if (p.nodes == new_path.nodes){
            add = false;
            break;
         }
      }
      if (add){
         P[i].push_front(new_path);
      }
   }
   */

   int cum5 = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now()-start5).count();

   // cout<<"times"<<endl;
   // cout<<cum1<<endl;
   // cout<<cum2<<endl;
   // cout<<cum3<<endl;
   // cout<<cum4<<endl;
   // cout<<cum5<<endl;

   for (int i = 0; i <=N.size(); i++)
   {
      for (int qu = 0; qu <= capacity; qu++)
      {
         P[i].insert(P[i].end(), P_q[qu][i].begin(), P_q[qu][i].end());
      }
   }

   return P;
}

list<SimpleRoute> GENROUTE(
   double z_ub,
   int Delta,
   double gamma,
   int h,
   bool &terminated,
   double &gamma_guarantee,
   VRP &vrp,
   DualSolution &sol,
   vector<vector<double>> penalized_distances
){
   int total_routes = 0;
   //int negative_routes = 0;
   int len_N = vrp.len_N();
   vector<int> quantities = vrp.quantities;
   vector<int> N = vrp.N;
   vector<vector<double>> &distance_dict = penalized_distances;
   vector<vector<double>> &geo_distance = vrp.geo_distance;
   int capacity = vrp.capacities[h - len_N];



   cout<<"     Generating paths"<<endl;

   double inf = numeric_limits<double>::infinity();
   double gamma_guarantee_1;
   GenPath P_l = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, "left", terminated, gamma_guarantee_1,inf);
   GenPath P_l_small = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, "left", terminated, gamma_guarantee_1,1);
   for (int i = 0; i<N.size(); i++){
      P_l[i].insert(P_l[i].end(), P_l_small[i].begin(), P_l_small[i].end());
   }

   // Debug
   if (h==5 && gamma == 5){
      //print_paths(P_l[1]);
   }
   double gamma_guarantee_2;
   GenPath P_r = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, "right", terminated, gamma_guarantee_2,inf);
   GenPath P_r_small = GENPATH(Delta, gamma, h, capacity, N, quantities, distance_dict, "right", terminated, gamma_guarantee_2,1);
   for (int i = 0; i<N.size(); i++){
      P_r[i].insert(P_r[i].end(), P_r_small[i].begin(), P_r_small[i].end());
   }


   cout<<"     Done with the paths"<<endl;
   cout<<"     Terminated: "<<terminated<<endl;
   if (h==5 && gamma == 5){
      //print_paths(P_r[1]);
      //p_v_v(distance_dict);
   }
   gamma_guarantee = min(gamma_guarantee_1,gamma_guarantee_2);
   // We order the paths by their cost
   auto cost_comparison = [](Path a, Path b) { return a.cost < b.cost; };
   for (int i = 0; i<len_N; i++){
      P_l[i].sort(cost_comparison);
      P_r[i].sort(cost_comparison);
   }
   cout<<"Done with the ordering"<<endl;
   /*
   if (Delta == std::numeric_limits<int>::max()){
      cout<<"len P_l: "<<P_l.size()<<endl;
      cout<<"len P_r: "<<P_l.size()<<endl;
   }
   */

   // for (int i = 0; i<len_N; i++){
   //    cout<<"Left Paths "<<i<<":"<<endl;
   //    print_paths(P_l[i]);
   //    cout<<"Right Paths "<<i<<":"<<endl;
   //    print_paths(P_r[i]);
   //
   // }
   //print_paths(P_l[len_N-1]);
   //print_paths(P_r[len_N-1]);

   vector<list<Route>> T(len_N, list<Route>(0));
   vector<list<Route>> R(len_N, list<Route>(0));

   // Set of added pairs of routes
   vector<set<std::tuple<int, int>>> added(len_N);

   // Initialize the routes
   for (int i=0; i<len_N; i++){
      added[i].insert(std::make_tuple(-1, -1));
      if (P_l[i].size()>0 && P_r[i].size()>0){
         Route init;
         init.path_l = P_l[i].begin();
         init.path_r = P_r[i].begin();
         init.cost = (init.path_l->cost)+ (init.path_r->cost);
         init.index_l = 0;
         init.index_r = 0;
         init.l_lb = (init.path_l->lower_bound);
         init.r_lb = (init.path_r->lower_bound);
         T[i].push_front(init);
         added[i].insert(std::make_tuple(0, 0));
      }
   }
   // Define infinity
   int iterations = 0;
   double cost_min;
   while(true){
      ++iterations;
      // Check the route with the smallest cost
      int arg_min = -1;
      cost_min = inf;
      for(int i = 0; i < len_N; i++){
         if (T[i].size() > 0){
            double cost = T[i].front().cost;
            if (cost<cost_min){
               cost_min = cost;
               arg_min = i;
            }
         }
      }
      // Break if no more routes
      if (arg_min == -1)
         break;
      // Break if exceed max cost
      if (cost_min>gamma){
         break;
      }
      // Check the first route and pop it
      Route r_star = T[arg_min].front();
      T[arg_min].pop_front();

      // Fill the load
      int load_l = r_star.path_l->load;
      int load_r = r_star.path_r->load;
      int total_load = load_l + load_r - quantities[arg_min];
      r_star.load = total_load;

      //Fill the median
      r_star.median = arg_min;
      bool valid = true;
      if (total_load > capacity){
         // Check if violate capacity
         valid = false;
      } else if (((double) load_l < (double)total_load/2.0) ||
                  ((double) load_r < (double)total_load/2.0) ||
                  ((double) load_l > (double)total_load/2.0 +(double) quantities[arg_min]) ||
                  ((double) load_r > (double)total_load/2.0 +(double) quantities[arg_min])){
         // Check if violate median
         valid = false;
      } else {
         // Check if violate empty intersection
         bit_set intersection = (r_star.path_l->nodes)&(r_star.path_r->nodes);
         bit_set comparison;
         comparison[arg_min] = 1;
         if (intersection != comparison){
            valid = false;
         }
      }
      if (valid){
         // Calculate the union of the nodes in the path
         r_star.nodes = (r_star.path_l->nodes)|(r_star.path_r->nodes);

         for (int i = 0; i<len_N; i++){
            for (auto r = R[i].begin(); r!=R[i].end(); ++r){
               if ((r->nodes) == (r_star.nodes)){
                  valid = false;
                  break;
               }
            }
            if (!valid)
               break;
         }
      }
      // Push if valid
      if (valid){
         R[arg_min].push_back(r_star);

         total_routes++;
      }
      // We add new routes
      vector<Route> new_routes;
      if (r_star.index_l + 1 < (int)P_l[arg_min].size()){
         // Check if not added already
         if (added[arg_min].count(std::make_tuple(r_star.index_l + 1, r_star.index_r)) == 0) {
            Route n_route;
            n_route.path_l = std::next(r_star.path_l,1);
            n_route.path_r = r_star.path_r;
            n_route.cost = (n_route.path_l->cost)+ (n_route.path_r->cost);
            n_route.l_lb = (n_route.path_l->lower_bound);
            n_route.r_lb = (n_route.path_r->lower_bound);
            n_route.index_l = r_star.index_l + 1;
            n_route.index_r = r_star.index_r;
            new_routes.push_back(n_route);
	    added[arg_min].insert(std::make_tuple(n_route.index_l, n_route.index_r));
         }
      }
      if (r_star.index_r + 1 < (int)P_r[arg_min].size()){
         if (added[arg_min].count(std::make_tuple(r_star.index_l, r_star.index_r + 1)) == 0){
            Route n_route;
            n_route.path_r = std::next(r_star.path_r,1);
            n_route.path_l = r_star.path_l;
            n_route.cost = (n_route.path_l->cost)+ (n_route.path_r->cost);
            n_route.l_lb = (n_route.path_l->lower_bound);
            n_route.r_lb = (n_route.path_r->lower_bound);
            n_route.index_r = r_star.index_r + 1;
            n_route.index_l = r_star.index_l;
            new_routes.push_back(n_route);
            added[arg_min].insert(std::make_tuple(n_route.index_l, n_route.index_r));
         }
      }
      // Order by cost
      std::sort (new_routes.begin(), new_routes.end(), compare_routes);
      // Insert them
      int loc = 0;
      auto p = T[arg_min].begin();
      while (p!=T[arg_min].end()){
         if (loc == (int) new_routes.size())
            break;
         if ((p->cost)>new_routes[loc].cost){
            T[arg_min].insert(p,new_routes[loc]);
            loc++;
         } else {
            ++p;
         }
      }
      // Finish inserting
      while (loc != (int) new_routes.size()){
         T[arg_min].insert(p,new_routes[loc]);
         loc++;
      }
   }


   cout<<"     Reconstructing SimpleRoutes"<<endl;
   // Take the routes and reconstruct them, since we do not want pointers to
   // iterators but the paths themselves
   std::list<SimpleRoute> SimpleRoutes;

   for (int i = 0; i < len_N; i++){
      for (auto route:R[i]){
         SimpleRoute copy;
         // Add first path
         for (auto it = (route.path_l->path).begin(); it!=(route.path_l->path).end(); it++){
            copy.path.push_back(*it);
         }
         // Add second path in opposite order (make sure not to add first element)
         for (auto it = (route.path_r->path).rbegin(); it!=(route.path_r->path).rend(); it++){
            if (it!=(route.path_r->path).rbegin())
               copy.path.push_back(*it);
         }
         copy.index_l = route.index_l;
         copy.index_r = route.index_r;
         copy.cost = route.cost;
         copy.l_lb = route.l_lb;
         copy.r_lb = route.r_lb;
         copy.load = route.load;
         copy.median = route.median;
         copy.truck = h;

         // Fill the geo_cost
         copy.geo_cost = 0;
         auto it1 = copy.path.begin();
         auto it2 = copy.path.begin();
         it2++;
         while(it2!=copy.path.end()){
            copy.geo_cost+=geo_distance[*it1][*it2];
            it1++;
            it2++;
         }
         if (copy.l_lb > copy.cost + pow(10,-13)*z_ub || copy.r_lb > copy.cost + pow(10,-13)*z_ub){
            print_sroute(copy);
            cout<<"LB not working"<<endl;
            cout<<copy.l_lb - copy.cost<<endl;
            cout<<copy.r_lb - copy.cost<<endl;
            throw "The lower bounds are not properly working";
         }
         SimpleRoutes.push_back(copy);
      }
   }

   cout<<"     Total of routes: "<<SimpleRoutes.size()<<endl;
   int total_neg = 0;
   for (auto r:SimpleRoutes){
      if (r.cost < -pow(10,-13) * z_ub){
         total_neg+=1;
         //print_sroute(r);
      }
   }
   cout<<"     Total of negative routes: "<<total_neg<<endl;

   // // We are also going to add all of the paths in the GENPATH by simply adding
   // // the house node at the end. This avoids convergence problems found in the code
   // for (int i = 0; i < len_N; i++){
   //    for (auto path:P_l[i]){
   //       bool add = true;
   //       SimpleRoute n_route;
   //       for (auto it = path.path.begin(); it!=path.path.end(); it++){
   //          n_route.path.push_back(*it);
   //       }
   //       n_route.path.push_back(h);
   //       n_route.cost = path.cost + distance_dict[path.end][h];
   //       n_route.index_r = -1;
   //       n_route.index_l = -1;
   //       n_route.l_lb = path.lower_bound;
   //       n_route.r_lb = -inf;
   //
   //       n_route.load = path.load;
   //       n_route.median = -1;
   //       n_route.truck = h;
   //       n_route.geo_cost = 0;
   //       auto it1 = n_route.path.begin();
   //       auto it2 = n_route.path.begin();
   //       it2++;
   //       while(it2!=n_route.path.end()){
   //          n_route.geo_cost+=geo_distance[*it1][*it2];
   //          it1++;
   //          it2++;
   //       }
   //       if (n_route.l_lb > n_route.cost + pow(10,-13) * z_ub|| n_route.r_lb > n_route.cost + pow(10,-13) * z_ub){
   //          print_sroute(n_route);
   //          throw "The lower bounds are not properly working";
   //       }
   //       if (n_route.cost < -pow(10,-13) * z_ub){
   //          // It is not necessary to check whether the route has been added or not, check later
   //          for (auto sroute:SimpleRoutes){
   //             if (sroute.path == n_route.path){
   //                //print_sroute(sroute);
   //                add = false;
   //                cout<<"Warning, added route that already exists"<<endl;
   //                //throw "Warning, added route that already exists";
   //                break;
   //             }
   //          }
   //          if (add){
   //             SimpleRoutes.push_back(n_route);
   //          }
   //       }
   //    }
   // }
   // for (int i = 0; i < len_N; i++){
   //    for (auto path:P_r[i]){
   //       bool add = true;
   //       SimpleRoute n_route;
   //       n_route.path.push_back(h);
   //       for (auto it = path.path.rbegin(); it!=path.path.rend(); it++){
   //          n_route.path.push_back(*it);
   //       }
   //       n_route.cost = path.cost + distance_dict[h][path.end];
   //       n_route.index_r = -1;
   //       n_route.index_l = -1;
   //       n_route.r_lb = path.lower_bound;
   //       n_route.l_lb = -inf;
   //       n_route.load = path.load;
   //       n_route.median = -1;
   //       n_route.truck = h;
   //       n_route.geo_cost = 0;
   //       auto it1 = n_route.path.begin();
   //       auto it2 = n_route.path.begin();
   //       it2++;
   //       while(it2!=n_route.path.end()){
   //          n_route.geo_cost+=geo_distance[*it1][*it2];
   //          it1++;
   //          it2++;
   //       }
   //       if (n_route.l_lb > n_route.cost + pow(10,-13)* z_ub || n_route.r_lb > n_route.cost + pow(10,-13)*z_ub){
   //          print_sroute(n_route);
   //          throw "The lower bounds are not properly working";
   //       }
   //       if (n_route.cost < -pow(10,-13) * z_ub){
   //          // It is not necessary to check whether the route has been added or not, check later
   //          for (auto sroute:SimpleRoutes){
   //             if (sroute.path == n_route.path){
   //                add = false;
   //                cout<<"Warning, added route that already exists"<<endl;
   //                //throw "Warning, added route that already exists";
   //                break;
   //             }
   //          }
   //          if (add){
   //             SimpleRoutes.push_back(n_route);
   //          }
   //       }
   //    }
   // }

   // Add the simple route of h,n,h if there are no routes for a node n
   for (int i = 0; i < len_N; i++){
      if (R[i].size() == 0){
         SimpleRoute n_route;
         n_route.path.push_back(h);
         n_route.path.push_back(i);
         n_route.path.push_back(h);
         n_route.cost = distance_dict[i][h]+ distance_dict[h][i];
         n_route.index_r = -1;
         n_route.index_l = -1;
         n_route.load = quantities[i];
         n_route.median = i;
         n_route.truck = h;
         n_route.l_lb = -inf;
         n_route.r_lb = -inf;

         // Fill the geo_cost
         n_route.geo_cost = 0;
         auto it1 = n_route.path.begin();
         auto it2 = n_route.path.begin();
         it2++;
         while(it2!=n_route.path.end()){
            n_route.geo_cost+=geo_distance[*it1][*it2];
            it1++;
            it2++;
         }
         // Before pushing the route, make sure it is feasible
         if (capacity>=quantities[i]){
            SimpleRoutes.push_back(n_route);
         }
      }
   }





   return SimpleRoutes;

}

// This function takes a set of routes with given reduced costs and finds
// the lower bound they induce, as well as updates the reduced costs and
// returns a new version of mu and lambda
LowerBound lower_bound_M2(
   vector<list<SimpleRoute>> Routes,
   VRP &vrp,
   vector<double> mu,
   vector<double> lamb
){

   vector<int>& H = vrp.H;
   vector<int>& capacities = vrp.capacities;
   vector<int>& N = vrp.N;
   vector<int>& quantities = vrp.quantities;
   vector<int>& n_trucks = vrp.n_trucks;
   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();

   // Infinity
   double inf = numeric_limits<double>::infinity();

   // Update the costs of the routes
   for (auto & truck_routes:Routes){
      for (auto & route:truck_routes){
         route.cost = route.geo_cost;
         route.cost -= mu[route.truck - len_N];
         for (auto node:route.path){
            if (node<len_N){
               route.cost -= lamb[node];
            }
         }
      }
   }
   /*
   for (auto & truck_routes:Routes){
      for (auto & route:truck_routes){
         cout<<"Cost "<<route.cost<<" Geo "<<route.geo_cost<<endl;
      }
   }
*/
   //The vectors of minima
   vector<vector<double>> b(len_N, vector<double>(len_H, inf));
   vector<vector<SimpleRoute>> b_routes(len_N, vector<SimpleRoute>(len_H));

   for (int h = 0; h < len_H; h++){
      for (auto route:Routes[h]){
         for (auto node:route.path){
            // Verify the node is a farmer
            if (node<len_N){
               double new_b = route.cost / ((double) route.load) * ((double)quantities[node]);
               if (new_b < b[node][h]){
                  b[node][h] = new_b;
                  b_routes[node][h] = route;
               }
            }
         }
      }

   }


   // The second vector of minima
   vector<double> b_min(len_N);
   vector<SimpleRoute> b_min_routes(len_N);
   for (int n = 0; n < len_N; n++){
      std::vector<double>::iterator h_it = std::min_element(b[n].begin(), b[n].end());
      int h_m = std::distance(b[n].begin(), h_it);
      b_min[n] = b[n][h_m];
      b_min_routes[n] = b_routes[n][h_m];
   }

   //Calculate number of visits to each node
   vector<vector<int>> visits(len_N,vector<int>(len_N,0));
   for (int n = 0; n < len_N; n++){
      for(auto it_route = b_min_routes[n].path.begin();  it_route!=b_min_routes[n].path.end(); ++it_route){
         if (*it_route < len_N){
            visits[n][*it_route] += 1;
         }
      }
   }


   // Construct theta
   vector<double> theta(len_N,0);
   for (int j = 0; j<len_N; j++){
      for (int i=0; i<len_N; i++){
         theta[j] += ((double) quantities[i])/((double) b_min_routes[i].load)*(double)visits[i][j];
      }
   }


   // Construct rho
   vector<double> rho(len_H,0);
   for (int n=0; n<len_N; n++){
      rho[b_min_routes[n].truck - len_N] += ((double) quantities[n])/((double) b_min_routes[n].load);
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

// Given a set of routes, calculate the best lower bounds that they generate
// We initialize at mu and lamb
void lower_bound_optimizer_M2(
   int iterations,
   double z_ub,
   double epsilon,
   VRP &vrp,
   DualSolution &solution){

   vector<double> lamb = solution.lamb;
   vector<double> mu = solution.v;
   vector<list<SimpleRoute>>& Routes = solution.routes;


   vector<int>& H = vrp.H;
   vector<int>& capacities = vrp.capacities;
   vector<int>& N = vrp.N;
   vector<int>& quantities = vrp.quantities;
   vector<vector<double>>& geo_distance = vrp.geo_distance;
   vector<int>& n_trucks = vrp.n_trucks;

   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();
   // Define infinity
   double infinity = numeric_limits<double>::infinity();
   // Define u
   vector<double> u(len_N);

   // Vectors to store the optimal values
   vector<double> lamb_opt(len_N);
   vector<double> u_opt(len_N);
   vector<double> v_opt(len_H);
   // Here we store the values of the iterations
   vector<double> values;
   double max_val = -infinity;
   vector<double> save_theta;
   vector<double> save_rho;
   double g_den_1 = 0;
   double g_den_2 = 0;
   double g = 0;
   for (int iteration = 0; iteration<iterations; iteration++){
      // We pass these routes to the algorithm that calculates the lower bound
      LowerBound lb = lower_bound_M2(Routes, vrp, mu, lamb);

      // Check if the lower bound that we get improves
      if (lb.z_lb > max_val){
         max_val = lb.z_lb;
         u_opt = lb.u;
         v_opt = mu;
         lamb_opt = lamb;
      }
      if (iteration % 10 == 0){
         cout<<"  Bound "<<iteration<<": "<<lb.z_lb<<endl;
      }

      //cout<<lamb<<endl;
      //cout<<mu<<endl;
      save_theta = lb.theta;
      save_rho = lb.rho;
      values.push_back(lb.z_lb);
      // We calculate g for the step of the algorithm
      g_den_1 = 0;
      g_den_2 = 0;
      for (int i = 0; i< len_N; i++){
         g_den_1 += pow(lb.theta[i]-1.0,2);
      }
      for (int i = 0; i< len_H; i++){
         g_den_2 += pow(lb.rho[i]-n_trucks[i],2);
      }
      g = (z_ub - lb.z_lb)/(g_den_1 + g_den_2);
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
            //cout<<"New epsilon "<<epsilon<<endl;
            values.clear();
         }
         int n_increases = 0;
         for (int i = (int)grad.size()-5; i<(int)grad.size(); i++ ){
            n_increases += (int)(grad[i]>0);
         }
         if (n_increases >= 5){
            epsilon = epsilon*1.2;
            //cout<<"New epsilon "<<epsilon<<endl;
            values.clear();
         }
      }
      // Check if we have reached a zero gradient
      double gradient_norm = 0;
      for (int i = 0; i < len_N; i++){
         gradient_norm += pow(lb.theta[i] - 1.0,2);
      }
      for (int i = 0; i < len_H; i++){
         gradient_norm += pow(lb.rho[i] - n_trucks[i],2);
      }

      if (gradient_norm < pow(10.0, -20)){
         cout<<"Reached zero gradient"<<endl;
         break;
      }
   }

   DualSolution new_bound;
   solution.z_lb = max_val;
   solution.u = u_opt;
   solution.v = v_opt;
   solution.lamb = lamb_opt;
   cout<<"New opt "<< max_val<<endl;
}

// Generates at most Delta routes of cost at most gamma
// if terminate[h-len_N] is true, then the routes are not calculated
// paste_bound is the bound to add the Routes, if a route has a bound below
// paste_bound, it is not added
// reinitialize determines whether routes are pasted at the end of initial_sol.routes or not
TerminatingCondition generateTruckMinRoutes(
   double z_ub,
   double Delta,
   double gamma,
   int h,
   DualSolution& initial_sol,
   VRP& vrp,
   double bound,
   bool paste
){
   // Define infinity
   double inf = numeric_limits<double>::infinity();

   //Extract solution
   double initial_gamma_guarantee = inf;

   //Create reduced distances
   if (initial_sol.reduced_distances.size() == 0){
      initial_sol.calc_reduced_distances(vrp.geo_distance);
   }


   //Calculate lengths
   int len_N = vrp.len_N();
   int len_H = vrp.len_H();


   // Create penalized distance matrix
   vector<vector<double>> penalized_distances;
   if (vrp.penalized){
      penalized_distances = penalized_matrix(
         initial_sol.reduced_distances,
         vrp.penalties[h - len_N],
         len_N,
         len_H,
         vrp.penalty_factor
      );
   } else {
      penalized_distances = initial_sol.reduced_distances;
   }

   // Initialize the parameters
   bool terminated_initial = true;
   double truck_guarantee;
   cout<<"Generating Routes for Truck: "<<h<<endl;

   list<SimpleRoute> Routes = GENROUTE(z_ub, Delta, gamma, h, terminated_initial, truck_guarantee, vrp, initial_sol, penalized_distances);
   cout<<"  Truck guarantee: "<<truck_guarantee<<endl;

   // Delete routes
   if (!paste){
      initial_sol.routes[h - len_N].clear();
   }
   int new_routes_count = 0;
   for (auto route:Routes){
      if (route.cost < bound){
         initial_sol.routes[h - len_N].push_back(route);
         ++new_routes_count;
      }
   }

   TerminatingCondition condition;
   condition.terminated = terminated_initial;
   condition.gamma_guarantee = truck_guarantee;
   condition.new_routes = new_routes_count;

   cout<<"  Routes Found: "<<new_routes_count<<endl;

   return condition;
}




DualSolution optimize_lower_bound_M2(
   int sub_iterations,
   double z_ub,
   int Delta,
   int Delta_zero,
   double gamma,
   double gamma_zero,
   double epsilon,
   VRP &vrp,
   DualSolution &initial_sol
){

   vector<double> mu = initial_sol.v;
   vector<double> lamb = initial_sol.lamb;
   vector<double> u = initial_sol.u;
   vector<list<SimpleRoute>> &initial_routes = initial_sol.routes;
   double &initial_gamma_guarantee = initial_sol.gamma_guarantee;

   vector<int>& H = vrp.H;
   vector<int>& capacities = vrp.capacities;
   vector<int>& N = vrp.N;
   vector<int>& quantities = vrp.quantities;
   vector<vector<double>>& geo_distance = vrp.geo_distance;
   vector<int>& n_trucks = vrp.n_trucks;

   //Calculate lengths
   int len_N = N.size();
   int len_H = H.size();


   // Calculate reduced costs
   if (initial_sol.reduced_distances.size() == 0){
      initial_sol.calc_reduced_distances(vrp.geo_distance);
   }

   // Define infinity
   double inf = numeric_limits<double>::infinity();

   // We start by generating routes for all of the trucks
   vector<TerminatingCondition> initial_terminating_conditions(len_H);

   for (auto h:H){
      initial_terminating_conditions[h-len_N] = generateTruckMinRoutes(
         z_ub,
         Delta,
         gamma,
         h,
         initial_sol,
         vrp,
         inf,
         false
      );
   }
   // // Now we add routes that come from a previous iteration of the algorithm
   // cout<<"Not Adding routes from previous M2 iteration "<<endl;
   // if (initial_routes.size()>0 && false){
   //    // Temporal place to store the routes
   //    vector<list<SimpleRoute>> unique_initial_routes(len_H);
   //    for (int i = 0; i< len_H; i++){
   //       for (auto old_route:initial_routes[i]){
   //          bool add = true;
   //          for (auto new_route:Routes[i]){
   //             if (new_route.geo_cost == old_route.geo_cost){
   //                if (new_route.path == old_route.path){
   //                   add = false;
   //                   break;
   //                }
   //             }
   //          }
   //          if (add){
   //             unique_initial_routes[i].push_back(old_route);
   //          }
   //       }
   //    }
   //    for (int i = 0; i< len_H; i++){
   //       for (auto old_route:unique_initial_routes[i]){
   //          Routes[i].push_back(old_route);
   //       }
   //    }
   // }
   //
   // initial_routes = Routes;

   // Here we store the new dual solution
   DualSolution ds = initial_sol;
   // Whether all routes are terminated, that is, the current mu,lamb are feasible
   bool all_terminated = false;

   while(true){

      // Reinitialize parameters
      // ds.u = initial_sol.u;
      // ds.lamb = initial_sol.lamb;
      // ds.v = initial_sol.v;

      cout<<"New iteration of gradient ascent"<<endl;
      // We pass these routes to the algorithm that calculates the lower bound
      // Check not using the cost
      // We use the same parameters of the previous solution

      lower_bound_optimizer_M2(sub_iterations, z_ub, epsilon, vrp, ds);
      ds.calc_reduced_distances(geo_distance);

      // Vectors indicate whether routes have been terminated or not
      vector<TerminatingCondition> terminating_conditions(len_H);

      int Delta_zero_current = Delta_zero;
      double total_new_routes_count = 0;

      cout<<"Calculating negative routes"<<endl;

      // This will run until either all trucks are optimal or negative routes have
      // been found
      while(true){
         for (auto h:H){
            // Verify truck has not terminated
            if (!terminating_conditions[h-len_N].terminated){
               terminating_conditions[h-len_N] = generateTruckMinRoutes(
                  z_ub,
                  Delta_zero_current,
                  gamma_zero,
                  h,
                  ds,
                  vrp,
                  - pow(10,-13) * z_ub,
                  true
               );
            }
         }

         total_new_routes_count = 0;
         for (auto h:H){
            total_new_routes_count += terminating_conditions[h-len_N].new_routes;
         }

         int total_non_terminated = 0;
         for (int i = 0; i < len_H; i++){
            if (!terminating_conditions[i].terminated){
               ++total_non_terminated;
            }
         }

         if (total_new_routes_count > 0){
            cout<<"Found negative routes: "<<total_new_routes_count<<endl;
            cout<<"Rerunning algorithm with more routes"<<endl;
            break;
         } else{
            if (total_non_terminated == 0){
               all_terminated = true;
               break;
            } else{
               cout<<"Need more routes to guarantee feasibility"<<endl;
               cout<<"Trucks without enough routes: "<<total_non_terminated<<endl;
               Delta_zero_current = Delta_zero_current * 2;
               cout<<"Trying new delta: "<<Delta_zero_current<<endl;
            }
         }

      }


      // cout<<"Routes per truck"<<endl;
      // for (int i = 0; i< len_H; i++){
      //    cout<<i+len_N<<":"<<Routes[i].size()<<endl;
      // }

      if (total_new_routes_count == 0){
         // We make sure that the algorithm has terminated because there are no
         // more negative routes
         if (all_terminated){
            break;
         } else {
            cout<<"Error, should not be seeing this message"<<endl;
         }

      }
   }
   cout<<"Reached feasibility!"<<endl;
   return ds;
}

vector<DualSolution> construct_lower_bound(
   int iterations_grad_m1,
   int iterations_grad_m2,
   int iterations_m2,
   double z_ub,
   int Delta,
   int Delta_zero,
   int Delta_final,
   double gamma,
   double gamma_zero,
   double gamma_final,
   double epsilon,
   VRP &vrp
){

   vector<int>& H = vrp.H;
   vector<int>& capacities = vrp.capacities;
   vector<int>& N = vrp.N;
   vector<int>& quantities = vrp.quantities;
   vector<vector<double>>& geo_distance = vrp.geo_distance;
   vector<int>& n_trucks = vrp.n_trucks;


   vector<DualSolution> solutions(1+iterations_m2);

   // Define infinity
   double inf = numeric_limits<double>::infinity();
   // Debugging first lower bounds
   DualSolution sol = lower_bound_optimizer_M1(iterations_grad_m1, z_ub, epsilon, vrp);
   DualSolution old_sol;
   cout<<"Bound Christofides:"<<sol.z_lb<<endl;
   sol.routes.clear();
   sol.initialize_routes(vrp.len_H());



   int len_N = N.size();
   int len_H = H.size();



   for (int iter_2 = 0; iter_2<iterations_m2; iter_2++){
      old_sol = sol;
      cout<<"Started Iteration of Bound 2 No. :"<<iter_2<<endl;
      sol = optimize_lower_bound_M2(iterations_grad_m2, z_ub, Delta, Delta_zero, gamma, gamma_zero, epsilon, vrp, old_sol);
      save_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+"_iter_"+to_string(iter_2)+".json", old_sol, vrp);


      // Debugging
      if (iter_2 == 0){
         for (int i = 0; i < len_H; i++){
            for (auto r:old_sol.routes[i]){
               if(r.cost<-pow(10,-13)*z_ub){
                  print_sroute(r);
                  cout<<"Chris"<<endl;
                  throw "Found route that violates Christofides";
               }
            }
         }
      }
      cout<<"Bound M2:"<<sol.z_lb<<endl;
      //Save the solutions
      solutions[iter_2] = old_sol;
   }
   cout<<"Finished Bounding"<<endl;

   // We start by generating routes for all of the trucks
   cout<<"Calculating Final Routes"<<endl;
   vector<TerminatingCondition> terminating_conditions(len_H);
   for (auto h:H){
      terminating_conditions[h-len_N] = generateTruckMinRoutes(
         z_ub,
         Delta_final,
         gamma_final,
         h,
         sol,
         vrp,
         - pow(10,-13) * z_ub,
         false
      );

   }
   //p_v(sol.u);
   //p_v(sol.v);
   cout<<"Terminated_Final Routes"<<endl;

   // Now we add routes that come from a previous iteration of the algorithm
   // We will not check whether the route exists or not to improve efficiency
   /*

   if (sol.routes.size()>0){
      for (int i = 0; i< len_H; i++){
         for (auto old_route:(sol.routes)[i]){
            bool add = true;
            for (auto new_route:FinalRoutes[i]){
               if (new_route.path == old_route.path){
                  add = false;
                  break;
               }
            }
            if (add){
               FinalRoutes[i].push_back(old_route);
            }
         }
      }
   }*/
   // cout<<"Not Joining Final Routes to Old Routes"<<endl;
   // // We will not add them
   // if (sol.routes.size()>0 && false){
   //    // Temporal place to store the routes
   //    vector<list<SimpleRoute>> unique_initial_routes(len_H);
   //    for (int i = 0; i< len_H; i++){
   //       for (auto old_route:sol.routes[i]){
   //          bool add = true;
   //          /*
   //          for (auto new_route:FinalRoutes[i]){
   //             if (new_route.geo_cost == old_route.geo_cost){
   //                if (new_route.path == old_route.path){
   //                   add = false;
   //                   break;
   //                }
   //             }
   //          }
   //          */
   //          if (add){
   //             unique_initial_routes[i].push_back(old_route);
   //          }
   //       }
   //    }
   //    for (int i = 0; i< len_H; i++){
   //       for (auto old_route:unique_initial_routes[i]){
   //          FinalRoutes[i].push_back(old_route);
   //       }
   //    }
   // }

   // cout<<"Routes per truck"<<endl;
   // for (int i = 0; i< len_H; i++){
   //    cout<<i+len_N<<":"<<FinalRoutes[i].size()<<endl;
   // }


   // Save the solution in the return argument
   solutions[iterations_m2] = sol;
   save_dual_solution(vrp.folder+"dual_solutions/"+vrp.name+"_iter_"+to_string(iterations_m2)+".json", sol, vrp);
   cout<<"Done"<<endl;
   return solutions;
}

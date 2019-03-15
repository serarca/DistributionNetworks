#pragma once

struct RouteParameters {
   int Delta;
   double gamma;
   double z_ub;
   // Limit for lower bounds, starts at 0
   int lb_limit;
   // Limit for the whole route, starts at 1
   int route_limit;
};

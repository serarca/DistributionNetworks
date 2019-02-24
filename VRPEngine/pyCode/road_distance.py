import urllib2
import json
import copy
# Calculates road distance between two (lat,lon) tuples, also returns time spent
# 
def road_distance(p0,p1, default_speed_mph = 6, detailed = False):
    if (not detailed):
        response_path = json.loads(urllib2.urlopen("http://127.0.0.1:5000/route/v1/driving/%f,%f;%f,%f"%(p0[1],p0[0],p1[1],p1[0])).read())
    else:
        response_path = json.loads(urllib2.urlopen("http://127.0.0.1:5000/route/v1/driving/%f,%f;%f,%f?steps=true"%(p0[1],p0[0],p1[1],p1[0])).read())

    #Calculate distance
    dist_road = response_path['routes'][0]['distance']
    response_closest_start = json.loads(urllib2.urlopen('http://127.0.0.1:5000/nearest/v1/foot/%f,%f'%(p0[1],p0[0])).read())
    dist_start = response_closest_start['waypoints'][0]['distance']
    response_closest_end = json.loads(urllib2.urlopen('http://127.0.0.1:5000/nearest/v1/foot/%f,%f'%(p1[1],p1[0])).read())
    dist_end = response_closest_end['waypoints'][0]['distance']
    total_dist = dist_road + dist_start + dist_end

    # Calculate time
    time_road = response_path['routes'][0]['duration']
    mile = 1.60934*1000
    hour = 60 * 60
    time_start = dist_start / mile / default_speed_mph * hour
    time_end = dist_end / mile / default_speed_mph * hour
    total_time = time_road + time_start + time_end

    result = {
        "distance":total_dist,
        "time":total_time
    }
    
    # Calculate the partition of the result
    if (detailed):
        detailed_distance = {"access":dist_start + dist_end}
        detailed_time = {"access":time_start + time_end}
        for leg in response_path['routes'][0]['legs']:
            for step in leg['steps']:
                surface = step["ref"]
                time = step["duration"]
                distance = step["distance"]
                if not (surface in detailed_distance.keys()):
                    detailed_distance[surface] = 0
                    detailed_time[surface] = 0
                detailed_distance[surface]+=distance
                detailed_time[surface]+=time
    
        result.update({
                "detailed_distance":detailed_distance,
                "detailed_time":detailed_time,
            })
        
    return result

def waypoints(p0,p1):
    response_path = json.loads(urllib2.urlopen("http://127.0.0.1:5000/route/v1/driving/%f,%f;%f,%f?steps=true&geometries=geojson&overview=full"%(p0[1],p0[0],p1[1],p1[0])).read())
    pol = response_path['routes'][0]['geometry']['coordinates']
    new_pol = copy.deepcopy(pol)
    for i in range(len(new_pol)):
        new_pol[i][0] = pol[i][1]
        new_pol[i][1] = pol[i][0]


    return new_pol


def path_distance(path):
    total_dist = 0
    for i in range(len(path) - 1):
        total_dist += road_distance(path[i],path[i+1])["distance"]

    return total_dist


import datetime
from dateutil import parser
from geopandas import GeoDataFrame
import pandas as pd
import numpy as np
from shapely.geometry import Point
import urllib
import json
import geopy.distance


# Calculates a distance matrix given a data frame
def geo_distance_matrix(df):
    df = df.reset_index(drop=True).copy()
    d = np.zeros([len(df),len(df)])
    for i,row_i in df.iterrows():
        for j,row_j in df.iterrows():
            if i!=j:
                d[i,j] = geopy.distance.vincenty(df.iloc[i].values,df.iloc[j].values).m
    return d


# Takes a data frame and calculates the distance (in m) and time (in min) between points
def gmaps_matrix(df, key):
    def format_distance(d):
        sp = d.split(" ")
        value = float(sp[0])
        if sp[1] == u'km':
            value = value * 1000
        return value

    def format_time(t):
        sp = t.split(" ")
        value = 0
        i = 0
        while i<len(sp):
            if sp[i+1]==u'min' or sp[i+1]==u'mins':
                value += float(sp[i])
            if sp[i+1]==u'hour' or sp[i+1]==u'hours':
                value += float(sp[i])*60
            i+=2
        return value
    indices = df.index.values
    sentence = np.char.array([str(a) for a in df.lat.values])+','+np.char.array([str(a) for a in df.lon.values])
    sentence = '|'.join(sentence)
    # Construct query
    query = ('https://maps.googleapis.com/maps/api/distancematrix/json?origins=' + sentence + 
             '&destinations=' + sentence + '&key=' + key)
    # Query the API
    f = urllib.urlopen(query)
    myfile = f.read()
    # Format output
    output = json.loads(myfile)
    n = len(output['rows'])
    distances = np.zeros([n,n])
    times = np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            distances[i,j] = format_distance(output['rows'][i]['elements'][j]['distance']['text'])
            times[i,j] = format_time(output['rows'][i]['elements'][j]['duration']['text'])
    df_distances = pd.DataFrame(data = distances, index = indices, columns = indices)
    df_times = pd.DataFrame(data = times, index = indices, columns = indices)
    return {'dist':df_distances,'time':df_times}
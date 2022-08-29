import geopy.distance
import numpy as np
from math import sin, cos, sqrt, atan2, radians

infile = 'geo_dis.txt'

def cal_dist(lat1, long1, lat2, long2):
    # approximate radius of earth in km
    R = 6373.0
    lat1 = radians(lat1)
    long1 = radians(long1)
    lat2 = radians(lat2)
    long2 = radians(long2)
    dlon = long2 - long1
    dlat = lat2 - lat1
    a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    distance = R * c
    return distance

HEADER = True
locs = []
names = []
with open(infile, 'r') as f:
    for l in f:
        if HEADER:
            HEADER = False
            continue
        items = l.split()
        lat = float(items[2])
        long = float(items[3])
        locs.append((lat, long))
        names.append(items[0])
print(locs)

matrix = np.zeros((len(locs), len(locs)))
for i in range(len(locs)):
    for j in range(len(locs)):
        if i == j:
            continue
        lat1, long1 = locs[i]
        lat2, long2 = locs[j]
        dis = geopy.distance.distance((lat1, long1), (lat2, long2)).km
        print(names[i] + ' ' + names[j] + ': ' + str(dis))     
        dis1 = cal_dist(lat1, long1, lat2, long2)
        print(names[i] + ' ' + names[j] + ': ' + str(dis1))
        matrix[i, j] = dis

print(names)
print(matrix)

np.savetxt('matrix.csv', matrix, delimiter='\t', header='\t'.join(names))

import math
from geographiclib.geodesic import Geodesic
from typing import Tuple, List

def get_bearing(pt1, pt2):
    bearing = Geodesic.WGS84.Inverse(pt1.lat, pt1.lon, pt2.lat, pt2.lon)['azi1']
    return normalize_bearing(bearing)

def normalize_bearing(bearing):
    if bearing < 0:
        return bearing + 360
    return bearing

def get_point_at_distance_and_bearing(point, distance, bearing):
    R = 6378100 #Radius of the Earth
    brng = math.radians(bearing) #Bearing is 90 degrees converted to radians.
    d = distance #Distance in km

    #lat2  52.20444 - the lat result I'm hoping for
    #lon2  0.36056 - the long result I'm hoping for.

    lat1 = math.radians(point.lat) #Current lat point converted to radians
    lon1 = math.radians(point.lon) #Current long point converted to radians

    lat2 = math.asin( math.sin(lat1)*math.cos(d/R) +
        math.cos(lat1)*math.sin(d/R)*math.cos(brng))

    lon2 = lon1 + math.atan2(math.sin(brng)*math.sin(d/R)*math.cos(lat1),
                math.cos(d/R)-math.sin(lat1)*math.sin(lat2))

    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)


    return Point(lat2, lon2)

def generate_triangle_points(center, radius):
    return [get_point_at_distance_and_bearing(center, radius, bearing) for bearing in [0,120,240]]


def great_circle_distance(ptA, ptB):
    # Degrees to radians
    phi1    = math.radians(ptA.lat)
    lambda1 = math.radians(ptA.lon)

    phi2    = math.radians(ptB.lat)
    lambda2 = math.radians(ptB.lon)

    delta_lambda = math.fabs(lambda2 - lambda1)

    central_angle = \
        math.atan2 \
        (
            # Numerator
            math.sqrt
            (
                # First
                math.pow
                (
                    math.cos(phi2) * math.sin(delta_lambda)
                    , 2.0
                )
                +
                # Second
                math.pow
                (
                    math.cos(phi1) * math.sin(phi2) -
                    math.sin(phi1) * math.cos(phi2) * math.cos(delta_lambda)
                    , 2.0
                )
            ),
            # Denominator
            (
                math.sin (phi1) * math.sin(phi2) +
                math.cos (phi1) * math.cos(phi2) * math.cos(delta_lambda)
            )
        )

    R = 6371009 # m
    return R * central_angle

class Point:
    def __init__(self, lat, lon):
        self.lat = lat
        self.lon = lon
    
    def to_tuple(self):
        return (self.lat, self.lon)

    def __repr__(self) -> str:
        return f'[{self.lat}, {self.lon}]'

    def __eq__(self, other):
        return self.lat == other.lat and self.lon == other.lon

class Beacon:
    def __init__(self, point, limits):
        self.point = point
        self.limits = limits

    def __eq__(self, other):
        return self.point == other.point and self.limits == other.limits

    def is_within_limits(self, point):
        dist = great_circle_distance(self.point, point)
        return dist >= self.limits[0] and dist <= self.limits[1]

    def get_limits(self):
        return [Limit(self, True), Limit(self, False)]

    def get_upper_limit(self):
        return Limit(self, False)

    def get_lower_limit(self):
        return Limit(self, True)

    def get_limits_with_distances_from_point(self, point):
        dist = great_circle_distance(self.point, point)
        return [(l, abs(dist - l.get_radius()))  for l in self.get_limits()]

    def get_nearest_limit_with_distance_from_point(self, point):
        return min(self.get_limits_with_distances_from_point(point), key=lambda l: l[1])

class Limit:
    def __init__(self, beacon, lower=True):
        self.beacon = beacon
        self.limit_index = 0 if lower else 1
    
    def get_radius(self):
        return self.beacon.limits[self.limit_index]

    def get_bearing_of_point(self, pt):
        return get_bearing(self.beacon.point, pt)

    def __eq__(self, other):
        return self.limit_index == other.limit_index and self.beacon == other.beacon

class IntersectionPoint:
    def __init__(self, pt, limit1, limit2):
        self.point = pt
        self.limit1 = limit1
        self.limit2 = limit2
    
    def __eq__(self, other):
        return self.point == other.point and self.limit1 == other.limit1 and self.limit2 == other.limit2

class BoundSquare:
    def __init__(self, center, dimension):
        self.north = get_point_at_distance_and_bearing(center, dimension/2, 0).lat
        self.south = get_point_at_distance_and_bearing(center, dimension/2, 180).lat
        self.west = get_point_at_distance_and_bearing(center, dimension/2, 270).lon
        self.east = get_point_at_distance_and_bearing(center, dimension/2, 90).lon


class Centroid:
    def __init__(self, point, radius):
        self.point = point
        self.radius = radius
        self.n = 0
    
    def new_empty(): 
        return Centroid(Point(0,0), 0)

    def new_from_points(points):
        c = Centroid.new_empty()
        c.add_points(points)
        c.divide()
        return c

    def add_points(self, points: List[Point]):
        for pt in points:
            self.add_point(pt)

    def add_point(self, point):
        self.point.lat += point.lat
        self.point.lon += point.lon
        self.n += 1

    def divide(self):
        self.point.lat /= self.n
        self.point.lon /= self.n
    
    def set_radius(self, points):
        self.radius = max([great_circle_distance(p, self.point) for p in points])

    def includes(self, pt):
        return great_circle_distance(pt, self.point) <= self.radius

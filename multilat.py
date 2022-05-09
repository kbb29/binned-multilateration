import math
from matplotlib.pyplot import magnitude_spectrum, plot
from scipy.optimize import minimize
from typing import Tuple, List
from plot import plotBeacons
from multilat_types import Point, Beacon, Limit, Centroid, IntersectionPoint, great_circle_distance, get_point_at_distance_and_bearing, generate_triangle_points
from limit_intersection import limit_intersection
from itertools import combinations, chain


def obfuscate_distance(dist, buckets):
    for bucket in buckets:
        if dist >= bucket[0] and dist < bucket[1]:
            return bucket

def simulate_service_distances(point, beacon_points, buckets):
    return [Beacon(bp, obfuscate_distance(great_circle_distance(bp, point), buckets)) for bp in beacon_points ]


def compute_intersection_points_for_beacon_pair(beacon_pair):
    (beacon1, beacon2) = beacon_pair
    ips  = limit_intersection( beacon1.get_lower_limit(), beacon2.get_upper_limit() )
    ips += limit_intersection( beacon1.get_lower_limit(), beacon2.get_lower_limit() )
    ips += limit_intersection( beacon1.get_upper_limit(), beacon2.get_upper_limit() )
    ips += limit_intersection( beacon1.get_upper_limit(), beacon2.get_lower_limit() )
    return ips

def compute_limit_intersections(beacons: List[Beacon]): # -> List[IntersectionPoint]:
    all_beacon_combinations = combinations(beacons, 2)
    return chain.from_iterable(map(compute_intersection_points_for_beacon_pair, all_beacon_combinations))

def is_limit_intersection_on_target(beacons: List[Beacon], limit_intersection: IntersectionPoint) -> bool:
    beacons_to_check = filter(lambda b: b != limit_intersection.limit1.beacon and b != limit_intersection.limit2.beacon, beacons)
    return all(b.is_within_limits(limit_intersection.point) for b in beacons_to_check)
    

def calculate_furthest_point(x, points):
    c = Centroid(Point(*x), 0)
    c.set_radius(points)
    return c.radius

def do_multilat(beacons: List[Beacon], tag=''):
    #compute all the points at which any of the limits of the beacons intersect
    all_intersections = compute_limit_intersections(beacons)

    #now filter the intersections to keep the ones which border the target area
    boundary_intersections = filter(lambda ip: is_limit_intersection_on_target(beacons, ip), all_intersections)
    bounds = list(map(lambda b: b.point,  boundary_intersections))
    
    #now calculate the centroid of these bounds
    bounds_centroid = Centroid.new_from_points(bounds)

    #now calculate the smallest circle which encloses the boundary points we mapped
    #this is an optimization based approach.
    #there are probably more efficient mathematical algorithms for doing this.
    result = minimize(calculate_furthest_point, bounds_centroid.point.to_tuple(), bounds, method='L-BFGS-B', options={'ftol': 1e-5, 'maxiter': 1e6})
    return Centroid(Point(*result.x), result.fun), bounds

if __name__ == '__main__':

    buckets = [[0,500], [500,1000], [1000, 2000], [2000, 5000], [5000, 10000], [10000, 20000]]
    center = Point(45,45)
    beacon_points = generate_triangle_points(center, 1100)
    for n,actual in enumerate([center, get_point_at_distance_and_bearing(center, 1300, 300)]): 
        beacons = simulate_service_distances(actual, beacon_points, buckets)
        centroid, bounds = do_multilat(beacons, tag=str(n))
            
        plotBeacons(beacons, actual=actual, preds=bounds, centroid=centroid, options={'tag': str(n)})
        
    beacon_points = generate_triangle_points(center, 1300) # + [get_point_at_distance_and_bearing(center, 15900, 180)]
    for n,actual in enumerate([Point(45.02245801238342, 45.03176042407004) ]):
        beacons = simulate_service_distances(actual, beacon_points, buckets)
        
        centroid, bounds = do_multilat(beacons, tag=str(n+4))
        
        plotBeacons(beacons, actual=actual, preds=bounds, centroid=centroid, options={'tag': str(n+4)})